# =============================================================
# Phase 1: Build the master cell-line index and lookup tables
# Produces THREE outputs:
#   1. master_cell_line_index.parquet   → ACH-ID + biological metadata
#   2. alias_lookup_table.parquet        → alias_name → ACH-ID (for HPA)
#   3. profile_lookup_table.parquet      → ProfileID → ACH-ID (for DepMap RNA)
# =============================================================

import pandas as pd

pd.options.mode.copy_on_write = True
import re

from src.utils import load_config, get_logger, get_path, ensure_processed_dir, save_parquet

log = get_logger(__name__)
cfg = load_config()


def load_sample_info() -> pd.DataFrame:
    """Load DepMap sample_info — the master cell-line registry."""
    path = get_path("paths.raw.nomenclature.depmap_sample")
    df = pd.read_csv(path)
    log.info(f"Loaded sample_info: {df.shape[0]:,} rows, {df.shape[1]} columns")
    return df


def load_cellosaurus_full() -> pd.DataFrame:
    """Load Cellosaurus WITHOUT filtering — filtering before the join was
    causing false negatives for legitimate cell lines (e.g. fibroblast lines
    like Hs-XXX.T series) that don't fall into 'Cancer'/'Transformed' categories."""
    path = get_path("paths.raw.nomenclature.cellosaurus")
    df = pd.read_csv(path)
    log.info(f"Loaded cellosaurus FULL (unfiltered): {df.shape[0]:,} rows")
    return df


def join_sample_to_cellosaurus(sample_df: pd.DataFrame, cello_df: pd.DataFrame) -> pd.DataFrame:
    merged = sample_df.merge(
        cello_df[["Accession (CVCL_xxxx)", "Synonyms", "Diseases", "Category",
                  "Species of origin", "Secondary accession number(s)"]],
        left_on="RRID",
        right_on="Accession (CVCL_xxxx)",
        how="left"
    )

    matched = merged["Accession (CVCL_xxxx)"].notna().sum()
    total = len(merged)
    log.info(f"Matched {matched:,}/{total:,} ({matched / total:.1%}) cell lines to Cellosaurus (unfiltered join)")

    target_categories = cfg["entity_resolution"]["target_categories"]
    merged = merged.copy()  # force independent object before assignment
    merged.loc[:, "is_target_category"] = merged["Category"].isin(target_categories)

    still_unmatched = merged[merged["Accession (CVCL_xxxx)"].isna()]
    if len(still_unmatched) > 0:
        log.warning(f"{len(still_unmatched):,} cell lines still unmatched after full-file join "
                    f"— likely engineered derivatives with no Cellosaurus entry")

    return merged


def build_alias_lookup(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Explode the Synonyms column into a long alias -> ACH-ID lookup table.
    This is what HPA's plain cell-line names (e.g. '143B') will match against.
    """
    strip_pattern = cfg["entity_resolution"]["strip_chars"]
    compiled_pattern = re.compile(strip_pattern)

    def make_stripped(name: str) -> str:
        return compiled_pattern.sub("", name).upper()

    records = []

    for _, row in merged_df.iterrows():
        ach_id = row["DepMap_ID"]

        primary_names = [
            row.get("cell_line_name"),
            row.get("stripped_cell_line_name"),
            row.get("CCLE_Name"),
        ]
        for name in primary_names:
            if pd.notna(name):
                clean_name = str(name).strip()
                records.append({
                    "alias": clean_name,
                    "DepMap_ID": ach_id,
                    "alias_stripped": make_stripped(clean_name),
                })

        synonyms = row.get("Synonyms")
        if pd.notna(synonyms):
            for alias in str(synonyms).split(";"):
                alias = alias.strip()
                if alias:
                    records.append({
                        "alias": alias,
                        "DepMap_ID": ach_id,
                        "alias_stripped": make_stripped(alias),
                    })

    # Single, one-shot construction — no post-hoc column assignment at all
    alias_df = pd.DataFrame.from_records(records).drop_duplicates(ignore_index=True)

    log.info(f"Built alias lookup: {len(alias_df):,} alias entries for "
             f"{alias_df['DepMap_ID'].nunique():,} unique cell lines")
    return alias_df


def build_profile_lookup() -> pd.DataFrame:
    """
    Load OmicsProfiles and filter to RNA profiles only.
    This maps ProfileID -> ACH-ID for joining the expression matrix.
    """
    path = get_path("paths.raw.nomenclature.depmap_omics")
    df = pd.read_csv(path)
    log.info(f"Loaded OmicsProfiles: {df.shape[0]:,} rows")

    rna_only = df[df["Datatype"] == "rna"].copy()
    log.info(f"Filtered to 'rna' datatype: {rna_only.shape[0]:,} profiles "
             f"covering {rna_only['ModelID'].nunique():,} unique cell lines")

    profile_lookup = rna_only[["ProfileID", "ModelID"]].rename(
        columns={"ModelID": "DepMap_ID"}
    )
    return profile_lookup


# Known derivative-naming patterns seen in the 22 unmatched rows
DERIVATIVE_PATTERNS = [
    r"^(.*?)(STAG2|TP53|KRAS)?(KO|NT)\d+$",  # A673STAG2KO16 -> A673
    r"^(.*?)-?ss\d+$",  # RPE1-ss48 -> RPE1
    r"^(.*?)_RPMI$",  # A375_RPMI -> A375
    r"^(.*?)GR$",  # HCC827GR -> HCC827
]


def extract_parent_name(name: str) -> str | None:
    """Attempt to strip known engineering/clone suffixes to recover the
    parent cell-line name, e.g. 'A673STAG2KO16' -> 'A673'."""
    for pattern in DERIVATIVE_PATTERNS:
        match = re.match(pattern, name, flags=re.IGNORECASE)
        if match and match.group(1):
            return match.group(1).strip("-_").upper()
    return None


def resolve_engineered_derivatives(unmatched_df: pd.DataFrame, alias_df: pd.DataFrame) -> pd.DataFrame:
    """
    For cell lines with no RRID and no Cellosaurus match, attempt to
    extract a parent line name and look it up in the alias table.
    Returns a dataframe: DepMap_ID -> resolved_parent_DepMap_ID (or None), method flag.
    """
    records = []

    for _, row in unmatched_df.iterrows():
        stripped_name = row["stripped_cell_line_name"]
        parent_name = extract_parent_name(stripped_name)

        resolved_id = None
        if parent_name:
            match = alias_df[alias_df["alias_stripped"] == parent_name.upper()]
            if not match.empty:
                resolved_id = match.iloc[0]["DepMap_ID"]

        records.append({
            "DepMap_ID": row["DepMap_ID"],
            "cell_line_name": row.get("cell_line_name"),
            "stripped_cell_line_name": stripped_name,
            "extracted_parent_name": parent_name,
            "resolved_parent_DepMap_ID": resolved_id,
            "resolution_method": "engineered_derivative" if resolved_id else "unresolved",
        })

    result_df = pd.DataFrame.from_records(records)

    resolved_count = result_df["resolved_parent_DepMap_ID"].notna().sum()
    log.info(f"Resolved {resolved_count}/{len(result_df)} engineered derivatives to a parent line")

    return result_df


def build_master_index():
    ensure_processed_dir()

    log.info("=" * 60)
    log.info("PHASE 1: Building master cell-line index")
    log.info("=" * 60)

    sample_df = load_sample_info()
    cello_full_df = load_cellosaurus_full()  # CHANGED: full, unfiltered

    merged_df = join_sample_to_cellosaurus(sample_df, cello_full_df)

    master_index = merged_df[[
        "DepMap_ID", "cell_line_name", "stripped_cell_line_name", "CCLE_Name",
        "RRID", "COSMICID", "primary_disease", "lineage", "lineage_subtype",
        "sample_collection_site", "primary_or_metastasis", "sex", "age",
        "default_growth_pattern", "model_manipulation",
        "Cellosaurus_NCIt_disease", "Diseases", "Category", "is_target_category"
    ]].copy()
    master_index = master_index.set_index("DepMap_ID")

    alias_lookup = build_alias_lookup(merged_df)
    profile_lookup = build_profile_lookup()

    # NEW: resolve remaining unmatched rows as engineered derivatives
    still_unmatched = merged_df[merged_df["Accession (CVCL_xxxx)"].isna()]
    if len(still_unmatched) > 0:
        derivative_resolution = resolve_engineered_derivatives(still_unmatched, alias_lookup)
        derivative_path = get_path("paths.processed.base") / "engineered_derivative_resolution.parquet"
        derivative_resolution.to_parquet(derivative_path, index=False)
        log.info(f"Saved derivative resolution table → engineered_derivative_resolution.parquet")

        truly_unresolved = derivative_resolution[derivative_resolution["resolution_method"] == "unresolved"]
        if len(truly_unresolved) > 0:
            log.info(
                f"{len(truly_unresolved)} cell lines have no Cellosaurus entry (e.g. lab-specific "
                f"subclones, media-adapted variants). This is expected — Cellosaurus does not catalogue "
                f"every individual derivative. These lines remain fully usable via their own DepMap ACH-ID; "
                f"they simply lack Cellosaurus-sourced synonym/disease metadata."
            )
            log.debug(f"Unresolved: {truly_unresolved['cell_line_name'].tolist()}")

    save_parquet(master_index, "paths.processed.master_index", log=log)
    save_parquet(alias_lookup, "paths.processed.alias_lookup", log=log)

    profile_lookup_path = get_path("paths.processed.base") / "profile_lookup_table.parquet"
    profile_lookup.to_parquet(profile_lookup_path, index=False)
    log.info(f"Saved {len(profile_lookup):,} rows → profile_lookup_table.parquet")

    log.success("Phase 1 complete: master index + alias lookup + profile lookup + derivative resolution ✓")
    return master_index, alias_lookup, profile_lookup


if __name__ == "__main__":
    build_master_index()

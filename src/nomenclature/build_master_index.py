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


def load_cellosaurus() -> pd.DataFrame:
    """Load Cellosaurus and filter to human cancer/transformed lines only."""
    path = get_path("paths.raw.nomenclature.cellosaurus")
    df = pd.read_csv(path)
    log.info(f"Loaded cellosaurus RAW: {df.shape[0]:,} rows")

    human_taxid = cfg["entity_resolution"]["human_taxid"]
    target_categories = cfg["entity_resolution"]["target_categories"]

    human_mask = df["Species of origin"].str.contains(human_taxid, na=False)
    category_mask = df["Category"].isin(target_categories)

    df_filtered = df[human_mask & category_mask].copy()
    log.info(f"Filtered to human + {target_categories}: {df_filtered.shape[0]:,} rows")
    return df_filtered

def join_sample_to_cellosaurus(sample_df: pd.DataFrame, cello_df: pd.DataFrame) -> pd.DataFrame:
    """Join on RRID (DepMap) = Accession (Cellosaurus)."""
    merged = sample_df.merge(
        cello_df[["Accession (CVCL_xxxx)", "Synonyms", "Diseases", "Category"]],
        left_on="RRID",
        right_on="Accession (CVCL_xxxx)",
        how="left"
    )

    matched = merged["Accession (CVCL_xxxx)"].notna().sum()
    total = len(merged)
    log.info(f"Matched {matched:,}/{total:,} ({matched/total:.1%}) cell lines to Cellosaurus")

    unmatched = merged[merged["Accession (CVCL_xxxx)"].isna()]
    if len(unmatched) > 0:
        log.warning(f"{len(unmatched):,} cell lines unmatched — will need fuzzy matching in a later step")

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

def build_master_index():
    """Main entry point — runs the full Phase 1 pipeline."""
    ensure_processed_dir()

    log.info("=" * 60)
    log.info("PHASE 1: Building master cell-line index")
    log.info("=" * 60)

    sample_df = load_sample_info()
    cello_df = load_cellosaurus()

    merged_df = join_sample_to_cellosaurus(sample_df, cello_df)

    # Keep the clean master index (metadata only, no exploded synonyms)
    master_index = merged_df[[
        "DepMap_ID", "cell_line_name", "stripped_cell_line_name", "CCLE_Name",
        "RRID", "COSMICID", "primary_disease", "lineage", "lineage_subtype",
        "sample_collection_site", "primary_or_metastasis", "sex", "age",
        "default_growth_pattern", "model_manipulation",
        "Cellosaurus_NCIt_disease", "Diseases", "Category"
    ]].copy()
    master_index = master_index.set_index("DepMap_ID")

    alias_lookup = build_alias_lookup(merged_df)
    profile_lookup = build_profile_lookup()

    # Save all three outputs
    save_parquet(master_index, "paths.processed.master_index", log=log)
    save_parquet(alias_lookup, "paths.processed.alias_lookup", log=log)

    profile_lookup_path = get_path("paths.processed.base") / "profile_lookup_table.parquet"
    profile_lookup.to_parquet(profile_lookup_path, index=False)
    log.info(f"Saved {len(profile_lookup):,} rows → profile_lookup_table.parquet")

    log.success("Phase 1 complete: 3 output files created ✓")
    return master_index, alias_lookup, profile_lookup


if __name__ == "__main__":
    build_master_index()

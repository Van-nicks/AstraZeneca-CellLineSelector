# =============================================================
# Phase 2, Step 3: Ingest HPA RNA expression (cell line panel)
# Bridges HPA "Cell line" names -> ACH-ID via:
#   hpa_desc (Cell line -> Cellosaurus ID)
#   depmap_sample.csv (RRID = Cellosaurus ID -> DepMap_ID)
# Chosen over parsing the 152K-row raw cellosaurus.csv Cross-references
# field, since depmap_sample.csv already has a direct RRID column
# covering 1,818/1,840 rows (98.8%) — simpler, reuses an existing file.
# Pivots long-format (Gene x Cell line x nTPM) to wide (ACH-ID x Gene).
# =============================================================

import pandas as pd

pd.options.mode.copy_on_write = True

from src.utils import load_config, get_logger, get_path, ensure_processed_dir, save_parquet

log = get_logger(__name__)
cfg = load_config()


def build_cellline_to_ach_map() -> pd.DataFrame:
    """
    Builds HPA 'Cell line' name -> ACH-ID map via two-hop bridge:
    hpa_desc (Cell line -> Cellosaurus ID) + depmap_sample (RRID -> ACH-ID).

    IMPORTANT: Rows with missing RRID/Cellosaurus ID are dropped from BOTH
    sides before merging. Pandas treats NaN == NaN as a match during merge,
    which caused a false fan-out (8 HPA cell lines x 22 NaN-RRID DepMap rows
    = 176 spurious rows) in initial testing. Confirmed via direct inspection
    that these are not real biological duplicates.
    """
    desc_path = get_path("paths.raw.nomenclature.hpa_desc")
    sample_path = get_path("paths.raw.nomenclature.depmap_sample")

    hpa_desc = pd.read_csv(desc_path, sep="\t")
    depmap_sample = pd.read_csv(sample_path)

    hpa_desc = hpa_desc.rename(columns={"Cellosaurus ID": "RRID"})

    hpa_desc_clean = hpa_desc[hpa_desc["RRID"].notna()]
    depmap_sample_clean = depmap_sample[depmap_sample["RRID"].notna()]

    dropped_hpa = len(hpa_desc) - len(hpa_desc_clean)
    if dropped_hpa > 0:
        log.info(f"Dropped {dropped_hpa} HPA cell lines with missing Cellosaurus ID "
                 f"before merge (avoids NaN-matches-NaN fan-out)")

    bridge = hpa_desc_clean[["Cell line", "RRID"]].merge(
        depmap_sample_clean[["DepMap_ID", "RRID"]], on="RRID", how="left"
    )

    matched = bridge["DepMap_ID"].notna().sum()
    log.info(f"HPA cell line -> ACH-ID bridge: {matched}/{len(bridge)} matched "
             f"({len(bridge) - matched} unmatched — dropped, likely non-DepMap cell lines)")

    dup_rrid_matches = bridge.groupby("Cell line").size()
    still_overmatched = dup_rrid_matches[dup_rrid_matches > 1]
    if len(still_overmatched) > 0:
        log.warning(f"{len(still_overmatched)} HPA cell lines still map to multiple ACH-IDs "
                    f"after NaN fix (real duplicate RRIDs, e.g. parent/derivative model pairs) — "
                    f"kept as-is; downstream pivot will produce one row per ACH-ID")

    return bridge.dropna(subset=["DepMap_ID"])[["Cell line", "DepMap_ID"]]


def ingest_hpa_expression() -> pd.DataFrame:
    rna_path = get_path("paths.raw.gene_expression.hpa_rna")

    log.info("Loading HPA RNA expression (long format, ~24M rows)...")
    hpa_rna = pd.read_csv(rna_path, sep="\t")
    log.info(f"Loaded: {hpa_rna.shape[0]:,} rows, {hpa_rna['Cell line'].nunique():,} cell lines, "
             f"{hpa_rna['Gene name'].nunique():,} genes")

    cellline_to_ach = build_cellline_to_ach_map()

    merged = hpa_rna.merge(cellline_to_ach, on="Cell line", how="inner")
    dropped_lines = hpa_rna["Cell line"].nunique() - merged["Cell line"].nunique()
    log.info(f"Retained {merged['Cell line'].nunique():,}/{hpa_rna['Cell line'].nunique():,} "
             f"cell lines after ACH-ID mapping ({dropped_lines} dropped)")

    log.info("Pivoting long format -> wide format (rows=ACH-ID, cols=Gene name, values=nTPM)...")
    wide = merged.pivot_table(index="DepMap_ID", columns="Gene name", values="nTPM", aggfunc="mean")

    log.info(f"Final HPA expression matrix: {wide.shape[0]:,} unique ACH-IDs x {wide.shape[1]:,} genes")
    return wide


def main():
    ensure_processed_dir()
    log.info("=" * 60)
    log.info("PHASE 2, STEP 3: Ingesting HPA RNA expression")
    log.info("=" * 60)

    hpa_df = ingest_hpa_expression()

    save_parquet(hpa_df, "paths.processed.hpa_expression", log=log)
    log.success("HPA RNA expression ingestion complete ✓")

    return hpa_df


if __name__ == "__main__":
    main()

# =============================================================
# Phase 2, Step 1: Ingest DepMap RNA expression matrix
# Re-indexes from ProfileID -> ACH-ID using profile_lookup_table.parquet
# Averages expression values for the 16 ACH-IDs with duplicate RNA profiles
# (documented design decision — no tiebreaker metadata exists in source data)
# =============================================================

import pandas as pd

pd.options.mode.copy_on_write = True

from src.utils import load_config, get_logger, get_path, ensure_processed_dir, save_parquet

log = get_logger(__name__)
cfg = load_config()


def load_profile_lookup() -> pd.DataFrame:
    path = get_path("paths.processed.base") / "profile_lookup_table.parquet"
    df = pd.read_parquet(path)
    log.info(f"Loaded profile_lookup_table: {len(df):,} rows")
    return df


def log_duplicate_profile_decision(profile_lookup: pd.DataFrame):
    """
    Documents the averaging decision for ACH-IDs with 2+ RNA profiles.
    Confirmed via inspection: OmicsProfiles.csv has no date/version/main-flag
    column to break ties — both profiles share identical ModelCondition.
    Design decision: average expression values (standard treatment for
    undocumented technical replicates).
    """
    dup_counts = profile_lookup.groupby("DepMap_ID")["ProfileID"].count()
    duplicated = dup_counts[dup_counts > 1]
    log.info(f"{len(duplicated)} ACH-IDs have multiple RNA profiles with no "
             f"tiebreaker metadata available — expression values will be averaged")


def ingest_depmap_expression() -> pd.DataFrame:
    expr_path = get_path("paths.raw.gene_expression.depmap")
    profile_lookup = load_profile_lookup()

    log.info("Loading DepMap expression matrix (large file — this may take a moment)...")
    expr_df = pd.read_csv(expr_path)
    expr_df = expr_df.rename(columns={"Unnamed: 0": "ProfileID"})
    log.info(f"Loaded expression matrix: {expr_df.shape[0]:,} profiles x {expr_df.shape[1] - 1:,} genes")

    log_duplicate_profile_decision(profile_lookup)

    merged = expr_df.merge(profile_lookup, on="ProfileID", how="inner")
    matched = merged["DepMap_ID"].notna().sum()
    log.info(f"Matched {matched:,}/{len(expr_df):,} profiles to ACH-ID")

    unmatched_profiles = set(expr_df["ProfileID"]) - set(merged["ProfileID"])
    if unmatched_profiles:
        log.warning(f"{len(unmatched_profiles)} expression profiles had no matching "
                    f"ACH-ID and were dropped")

    merged = merged.drop(columns=["ProfileID"])
    expression_by_ach = merged.groupby("DepMap_ID").mean(numeric_only=True)

    log.info(f"Final expression matrix: {expression_by_ach.shape[0]:,} unique ACH-IDs "
             f"x {expression_by_ach.shape[1]:,} genes (duplicates averaged)")

    return expression_by_ach


def main():
    ensure_processed_dir()
    log.info("=" * 60)
    log.info("PHASE 2, STEP 1: Ingesting DepMap RNA expression")
    log.info("=" * 60)

    expression_df = ingest_depmap_expression()

    save_parquet(expression_df, "paths.processed.depmap_expression", log=log)
    log.success("DepMap RNA expression ingestion complete ✓")

    return expression_df


if __name__ == "__main__":
    main()

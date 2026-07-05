# =============================================================
# Phase 2, Step 2: Ingest CCLE/Gygi proteomics matrix
# Already indexed directly by ACH-ID — no ProfileID mapping needed.
# Columns are KEPT as "UNIPROT (GENE)" (original format) because
# stripping to bare gene symbols caused duplicate column names —
# multiple distinct UniProt IDs (different proteins/isoforms) can
# share the same gene symbol (e.g. two entries both map to 'TPM1').
# A separate gene-symbol lookup (uniprot_to_gene map) is provided
# for cross-referencing with other omics layers, without collapsing
# or discarding any protein-level columns.
# =============================================================

import pandas as pd

pd.options.mode.copy_on_write = True

from src.utils import load_config, get_logger, get_path, ensure_processed_dir, save_parquet

log = get_logger(__name__)
cfg = load_config()


def build_uniprot_to_gene_map(columns) -> dict:
    """
    Builds a UniProtID -> GeneSymbol lookup from columns like
    'A0AV96 (RBM47)', WITHOUT renaming/collapsing the actual dataframe
    columns. Useful later for cross-referencing gene symbols across
    omics layers while keeping protein-level granularity intact here.
    """
    mapping = {}
    for col in columns:
        if "(" in col and col.endswith(")"):
            uniprot_id = col.split("(")[0].strip()
            gene_symbol = col.split("(")[-1].rstrip(")").strip()
            mapping[uniprot_id] = gene_symbol
    return mapping


def ingest_proteomics() -> pd.DataFrame:
    path = get_path("paths.raw.gene_expression.proteomics")

    log.info("Loading CCLE/Gygi proteomics matrix...")
    df = pd.read_csv(path)
    df = df.rename(columns={"Unnamed: 0": "DepMap_ID"})
    log.info(f"Loaded proteomics matrix: {df.shape[0]:,} cell lines x {df.shape[1] - 1:,} proteins")

    protein_cols = [c for c in df.columns if c != "DepMap_ID"]
    uniprot_to_gene = build_uniprot_to_gene_map(protein_cols)
    gene_counts = pd.Series(list(uniprot_to_gene.values())).value_counts()
    shared_genes = gene_counts[gene_counts > 1]
    if len(shared_genes) > 0:
        log.info(f"{len(shared_genes)} gene symbols map to multiple distinct UniProt "
                 f"entries (isoforms/paralogs) — columns kept as 'UNIPROT (GENE)' "
                 f"to avoid collisions")

    df = df.set_index("DepMap_ID")

    missing_pct = df.isna().mean().mean() * 100
    log.info(f"Average missingness across proteins: {missing_pct:.1f}%")

    dup_ach_ids = df.index[df.index.duplicated()].unique()
    if len(dup_ach_ids) > 0:
        log.warning(f"{len(dup_ach_ids)} duplicate ACH-IDs found — averaging")
        df = df.groupby(df.index).mean(numeric_only=True)
    else:
        log.info("No duplicate ACH-IDs — each cell line has exactly one profile")

    log.info(f"Final proteomics matrix: {df.shape[0]:,} unique ACH-IDs x {df.shape[1]:,} proteins")
    return df


def main():
    ensure_processed_dir()
    log.info("=" * 60)
    log.info("PHASE 2, STEP 2: Ingesting CCLE/Gygi proteomics")
    log.info("=" * 60)

    proteomics_df = ingest_proteomics()

    save_parquet(proteomics_df, "paths.processed.proteomics", log=log)
    log.success("Proteomics ingestion complete ✓")

    return proteomics_df


if __name__ == "__main__":
    main()

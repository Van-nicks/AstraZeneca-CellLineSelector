# AstraZeneca-CellLineSelector
Building a system that helps scientists find the "perfect" cell line for their specific experiment. If a scientist wants to target a specific mutation with a specific metabolic behaviour, this tool would help them find the cell lines that match those requirements.

## Progress: 5 / 14 raw files connected

### Phase 1 — Nomenclature Join (Complete)
Built a master `profile_lookup_table.parquet` mapping ProfileID -> ModelCondition
-> ACH-ID (DepMap_ID) across all cell lines, using:
- `depmap_sample` (9_DepMap_sample_info.csv)
- `cellosaurus` (7_cellosaurus.csv)
- `depmap_omics` (8_DepMap_OmicsProfiles.csv)

**Key design decision:** 16 ACH-IDs have multiple RNA profiles under identical
`ModelCondition` with no tiebreaker metadata in the source data (no date,
version, or "is main" flag exists). Decision: treat as technical replicates
and average their values downstream rather than arbitrarily discarding one.

Output: `data/processed/profile_lookup_table.parquet` (1,495 profile rows,
1,479 unique ACH-IDs)

### Phase 2 — Gene Expression Ingestion (In Progress)

#### Step 1: DepMap RNA Expression ✅
- Script: `src/ingestion/ingest_depmap.py`
- Source: `2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv`
- Re-indexed from ProfileID -> ACH-ID via profile_lookup_table
- Duplicate RNA profiles (16 ACH-IDs) averaged per Phase 1 decision
- Output: `data/processed/depmap_expression.parquet`
  — **1,479 rows x 53,961 genes**

#### Step 2: CCLE/Gygi Proteomics ✅
- Script: `src/ingestion/ingest_proteomics.py`
- Source: `4_Harmonized_MS_CCLE_Gygi_subsetted.csv`
- Already indexed directly by ACH-ID — no ProfileID mapping needed
- **Key design decision:** protein columns kept in original
  `UNIPROT (GENE)` format rather than renamed to bare gene symbols.
  396 gene symbols map to multiple distinct UniProt entries
  (isoforms/paralogs); renaming to gene-only caused column name
  collisions and a `pyarrow` write failure. A separate
  `uniprot_to_gene` lookup dict is built in-script for future
  cross-referencing without collapsing protein-level columns.
- Missingness: ~27.8% average (expected for mass-spec proteomics —
  not every protein detected in every cell line; no imputation
  applied yet)
- Output: `data/processed/proteomics.parquet`
  — **375 rows x 12,558 proteins**

### Remaining Files (Not Yet Ingested)
| File | Category | Notes |
|---|---|---|
| `hpa_rna` (1_4_hpa_rna_celline.tsv) | Gene expression | Likely indexed by cell line name, not ACH-ID — needs alias bridge via `hpa_desc` or `cellosaurus` |
| `geo` (3_GEOexpression.txt) | Gene expression | Not yet inspected |
| `hpa_desc` (11_hpa_rna_celline_description.tsv) | Nomenclature | Needed as bridge for HPA RNA |
| `geo_info` (10_GEOInfo.txt) | Nomenclature | Not yet inspected |
| `fusions` (5_OmicsFusionFilteredSupplementary.csv) | Gene properties | Not yet inspected |
| `mutations` (6_OmicsSomaticMutationsProfile.csv) | Gene properties | Not yet inspected |
| `metabolomics` (12_CCLE_metabolomics_20190502.csv) | Non-gene expression | Not yet inspected |
| `mirna` (13_CCLE_miRNA_20181103.gct) | Non-gene expression | Not yet inspected |
| `signatures` (14_OmicsGlobalSignatures.csv) | Non-gene expression | Not yet inspected |

## Data Versioning
Processed `.parquet` files are intentionally excluded from Git
(`.gitignore`) since they are large binary artifacts fully
regenerable from raw data + ingestion scripts. Raw data files are
also expected to be present locally (not committed). To reproduce
all processed outputs, run each ingestion script in
`src/ingestion/` in order, or a future `run_all.py` orchestrator.

## Config
All raw/processed file paths are defined centrally in `settings.yaml`
under `paths.raw.*` and `paths.processed.*`, keyed by category
(`nomenclature`, `gene_expression`, `gene_properties`,
`non_gene_expression`).
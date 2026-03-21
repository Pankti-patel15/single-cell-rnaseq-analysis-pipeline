# Single-cell RNA-seq Analysis Pipeline

This repository contains a Python pipeline for single-cell RNA-seq analysis using the public 10x Genomics PBMC 3k dataset.

The workflow starts from the downloaded 10x count matrix, performs standard preprocessing and clustering, identifies cluster-level marker genes, assigns broad immune cell type labels heuristically, and writes figures and summary tables for downstream interpretation.

## Dataset

The analysis uses the public PBMC 3k dataset from 10x Genomics:

- [PBMC 3k filtered gene-barcode matrices](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)

This was chosen because it is a real and widely used single-cell dataset that is small enough to run locally while still capturing the core steps of a practical scRNA-seq workflow.

## Workflow

The pipeline performs the following steps:

1. download and extract the 10x matrix files
2. calculate QC metrics
3. filter low-quality cells and low-detection genes
4. remove cells with elevated mitochondrial content
5. normalize counts and log-transform expression
6. select highly variable genes
7. run PCA, construct a nearest-neighbor graph, and generate a UMAP embedding
8. cluster cells with Leiden clustering
9. identify marker genes for each cluster
10. assign broad immune cell type labels based on canonical marker panels

## Results

In the current starter run, the pipeline produced:

- 2,700 cells before QC
- 2,698 cells after QC
- 13,714 genes retained after filtering
- 8 Leiden clusters

The resulting clusters were consistent with expected PBMC populations, including:

- T cells
- CD14 monocytes
- FCGR3A monocytes
- B cells
- NK cells
- platelets

Representative marker genes in the output included:

- `CD3D` for T-cell enriched clusters
- `S100A8`, `S100A9`, and `LYZ` for monocyte-enriched clusters
- `CD79A`, `CD74`, and `MS4A1` for B-cell enriched clusters
- `NKG7`, `PRF1`, and `GNLY` for NK-cell enriched clusters
- `PPBP` and `PF4` for platelets

## Output Files

After a successful run, the pipeline writes:

- `output/processed/cluster_summary.tsv`
- `output/processed/marker_genes.tsv`
- `output/reports/analysis_summary.md`
- `output/figures/qc_histograms.png`
- `output/figures/umap_clusters.png`
- `output/figures/umap_cell_types.png`
- `output/figures/marker_heatmap.png`

The full per-cell metadata and `.h5ad` object are generated locally but are typically not committed to Git because they are larger intermediate outputs.

## Project Structure

```text
.
├── run_pipeline.py
├── requirements.txt
├── pyproject.toml
├── src/single_cell_pipeline/
│   ├── analysis.py
│   ├── data.py
│   └── pipeline.py
└── tests/
```

## Running the Pipeline

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 run_pipeline.py
```

## Implementation Notes

This project is intended to show practical familiarity with standard single-cell RNA-seq analysis steps in Python rather than to introduce a novel method. The annotation step is heuristic and based on canonical immune markers, so the assigned labels should be interpreted as broad cell-type calls rather than final expert annotations.

## Testing

```bash
python3 -m unittest discover -s tests -v
```

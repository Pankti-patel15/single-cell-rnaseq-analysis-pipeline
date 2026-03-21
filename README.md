# Single-cell RNA-seq Analysis Pipeline

This repository contains an end-to-end Python pipeline for single-cell RNA-seq analysis on a real public dataset.

The workflow uses the 10x Genomics PBMC 3k dataset, starting from the raw 10x count matrix and proceeding through quality control, filtering, normalization, dimensionality reduction, clustering, marker gene analysis, and broad cell type annotation. The pipeline generates figures, summary tables, and a short report automatically.

## Dataset

The analysis uses the public PBMC 3k dataset from 10x Genomics:

- [PBMC 3k filtered gene-barcode matrices](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)

This dataset is widely used as a reference example for single-cell RNA-seq workflows and is well suited for demonstrating core preprocessing, clustering, and interpretation steps in a reproducible local analysis pipeline.

## Workflow

The pipeline performs the following steps:

1. downloads and extracts the 10x matrix files
2. calculates standard cell-level QC metrics
3. filters low-quality cells and low-detection genes
4. removes cells with elevated mitochondrial transcript content
5. normalizes counts and log-transforms expression
6. selects highly variable genes
7. runs PCA and constructs a nearest-neighbor graph
8. generates a UMAP embedding
9. performs Leiden clustering
10. identifies cluster-level marker genes
11. assigns broad immune cell type labels using canonical marker panels

## Results

In the current run, the pipeline produced:

- 2,700 cells before QC
- 2,698 cells after QC
- 32,738 genes before filtering
- 13,714 genes after filtering
- 8 Leiden clusters

The resulting clusters were consistent with expected PBMC populations, including:

- T cells
- CD14 monocytes
- FCGR3A monocytes
- B cells
- NK cells
- platelets

Representative marker genes recovered by the pipeline included:

- `CD3D` in T-cell enriched clusters
- `S100A8`, `S100A9`, and `LYZ` in monocyte-enriched clusters
- `CD79A`, `CD74`, and `MS4A1` in B-cell enriched clusters
- `NKG7`, `PRF1`, and `GNLY` in NK-cell enriched clusters
- `PPBP` and `PF4` in platelets

These outputs support the biological plausibility of the clustering and annotation steps.

## Output Files

After a successful run, the pipeline writes:

- `output/processed/cluster_summary.tsv`
- `output/processed/marker_genes.tsv`
- `output/reports/analysis_summary.md`
- `output/figures/qc_histograms.png`
- `output/figures/umap_clusters.png`
- `output/figures/umap_cell_types.png`
- `output/figures/marker_heatmap.png`

The full per-cell metadata table and `.h5ad` object are also generated locally for downstream analysis, but are typically not committed to Git because they are larger intermediate files.

## Repository Structure

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

This project is intended to demonstrate practical familiarity with standard single-cell RNA-seq analysis steps in Python. The cluster labels are heuristic and based on canonical immune markers, so they should be interpreted as broad cell type assignments rather than definitive annotations.

## Testing

```bash
python3 -m unittest discover -s tests -v
```

## License

MIT License. See `LICENSE`.

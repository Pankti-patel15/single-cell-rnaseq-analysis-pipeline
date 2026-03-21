from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .analysis import (
    CELL_TYPE_MARKERS,
    assign_cluster_cell_types,
    percent_mt_label,
    select_heatmap_markers,
    summarize_clusters,
    top_marker_table,
    write_summary_report,
)
from .data import download_dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="End-to-end single-cell RNA-seq analysis pipeline using a real public dataset."
    )
    parser.add_argument("--output-dir", default="output")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--max-mt-pct", type=float, default=15.0)
    parser.add_argument("--n-hvgs", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--resolution", type=float, default=0.5)
    return parser.parse_args()


def prepare_directories(output_dir: Path) -> dict[str, Path]:
    data_dir = output_dir / "data"
    processed_dir = output_dir / "processed"
    figures_dir = output_dir / "figures"
    reports_dir = output_dir / "reports"
    mpl_dir = output_dir / ".mplconfig"

    for path in (data_dir, processed_dir, figures_dir, reports_dir, mpl_dir):
        path.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    return {
        "data_dir": data_dir,
        "processed_dir": processed_dir,
        "figures_dir": figures_dir,
        "reports_dir": reports_dir,
    }


def require_scanpy():
    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "This single-cell project requires scanpy and related dependencies. "
            "Install them with: pip install -r requirements.txt"
        ) from exc
    return sc


def run_scanpy_workflow(args: argparse.Namespace, paths: dict[str, Path]):
    sc = require_scanpy()

    matrix_dir = download_dataset(paths["data_dir"])
    print(f"Loading 10x matrix from {matrix_dir}...")
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()

    cells_before, genes_before = adata.n_obs, adata.n_vars
    print(f"Loaded {cells_before} cells and {genes_before} genes before QC.")

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mt_pct].copy()
    cells_after, genes_after = adata.n_obs, adata.n_vars
    print(f"Retained {cells_after} cells and {genes_after} genes after QC.")

    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    for cell_type, markers in CELL_TYPE_MARKERS.items():
        valid = [marker for marker in markers if marker in adata.raw.var_names]
        if valid:
            sc.tl.score_genes(adata, gene_list=valid, score_name=f"{cell_type}_score", use_raw=True)

    annotated_obs = assign_cluster_cell_types(adata.obs)
    adata.obs = annotated_obs
    cluster_summary = summarize_clusters(adata.obs)
    marker_table = top_marker_table(adata)

    adata.write(paths["processed_dir"] / "processed_adata.h5ad")
    adata.obs.to_csv(paths["processed_dir"] / "cell_metadata.tsv", sep="\t")
    cluster_summary.to_csv(paths["processed_dir"] / "cluster_summary.tsv", sep="\t", index=False)
    marker_table.to_csv(paths["processed_dir"] / "marker_genes.tsv", sep="\t", index=False)

    return adata, cluster_summary, marker_table, cells_before, cells_after, genes_before, genes_after


def plot_qc(adata, destination: Path) -> None:
    mt_column = percent_mt_label(adata.obs)
    plot_df = adata.obs[["n_genes_by_counts", "total_counts", mt_column]].copy()
    plot_df.columns = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sns.histplot(plot_df["n_genes_by_counts"], bins=40, ax=axes[0], color="#1d3557")
    axes[0].set_title("Genes per Cell")
    sns.histplot(plot_df["total_counts"], bins=40, ax=axes[1], color="#457b9d")
    axes[1].set_title("UMI Counts per Cell")
    sns.histplot(plot_df["pct_counts_mt"], bins=40, ax=axes[2], color="#e76f51")
    axes[2].set_title("Mitochondrial Percentage")
    for ax in axes:
        ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(destination, dpi=300)
    plt.close()


def plot_umap(adata, color_column: str, destination: Path, title: str) -> None:
    plot_df = pd.DataFrame(
        {
            "UMAP1": adata.obsm["X_umap"][:, 0],
            "UMAP2": adata.obsm["X_umap"][:, 1],
            color_column: adata.obs[color_column].astype(str).values,
        }
    )

    plt.figure(figsize=(9, 7))
    sns.scatterplot(
        data=plot_df,
        x="UMAP1",
        y="UMAP2",
        hue=color_column,
        s=16,
        linewidth=0,
        palette="tab10",
    )
    plt.title(title)
    plt.tight_layout()
    plt.savefig(destination, dpi=300)
    plt.close()


def plot_marker_heatmap(adata, marker_table: pd.DataFrame, destination: Path) -> None:
    genes = [gene for gene in select_heatmap_markers(marker_table) if gene in adata.raw.var_names]
    if not genes:
        return

    sc = require_scanpy()
    plot_obj = sc.pl.heatmap(
        adata,
        var_names=genes,
        groupby="cluster_cell_type",
        use_raw=True,
        show=False,
        cmap="viridis",
    )
    fig = plot_obj["heatmap_ax"].figure
    fig.suptitle("Top Marker Genes Across Annotated Cell Types", y=1.02)
    fig.tight_layout()
    fig.savefig(destination, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    paths = prepare_directories(output_dir)

    print("Running single-cell RNA-seq pipeline...")
    adata, cluster_summary, marker_table, cells_before, cells_after, genes_before, genes_after = run_scanpy_workflow(
        args,
        paths,
    )

    print("Generating plots and summary report...")
    plot_qc(adata, paths["figures_dir"] / "qc_histograms.png")
    plot_umap(adata, "leiden", paths["figures_dir"] / "umap_clusters.png", "UMAP by Leiden Cluster")
    plot_umap(
        adata,
        "cluster_cell_type",
        paths["figures_dir"] / "umap_cell_types.png",
        "UMAP by Heuristic Cell Type Annotation",
    )
    plot_marker_heatmap(adata, marker_table, paths["figures_dir"] / "marker_heatmap.png")

    write_summary_report(
        paths["reports_dir"] / "analysis_summary.md",
        dataset_name="10x Genomics PBMC 3k",
        cells_before=cells_before,
        cells_after=cells_after,
        genes_before=genes_before,
        genes_after=genes_after,
        cluster_summary=cluster_summary,
        marker_table=marker_table,
    )

    print("\nPipeline completed successfully.")
    print(f"Outputs written to: {output_dir}")
    print(f"- Figures: {paths['figures_dir']}")
    print(f"- Results: {paths['processed_dir']}")
    print(f"- Report: {paths['reports_dir'] / 'analysis_summary.md'}")


if __name__ == "__main__":
    main()


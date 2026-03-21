from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


CELL_TYPE_MARKERS: dict[str, list[str]] = {
    "T_cells": ["CD3D", "CD3E", "IL7R", "LTB"],
    "NK_cells": ["NKG7", "GNLY", "PRF1", "GZMB"],
    "B_cells": ["MS4A1", "CD79A", "CD74", "HLA-DRA"],
    "CD14_monocytes": ["LYZ", "S100A8", "S100A9", "CTSS"],
    "FCGR3A_monocytes": ["FCGR3A", "LST1", "IFITM3", "TYMP"],
    "Dendritic_cells": ["FCER1A", "CST3", "HLA-DRA", "CLEC10A"],
    "Platelets": ["PPBP", "PF4", "SDPR", "NRGN"],
}


def assign_cluster_cell_types(
    obs: pd.DataFrame,
    cluster_column: str = "leiden",
    score_suffix: str = "_score",
) -> pd.DataFrame:
    score_columns = [col for col in obs.columns if col.endswith(score_suffix)]
    if not score_columns:
        raise ValueError("No marker score columns found for annotation.")

    per_cell_calls = obs[score_columns].idxmax(axis=1).str.replace(score_suffix, "", regex=False)
    annotated = obs.copy()
    annotated["predicted_cell_type"] = per_cell_calls

    cluster_map = {}
    for cluster, frame in annotated.groupby(cluster_column):
        winner = Counter(frame["predicted_cell_type"]).most_common(1)[0][0]
        cluster_map[cluster] = winner

    annotated["cluster_cell_type"] = annotated[cluster_column].map(cluster_map)
    return annotated


def summarize_clusters(obs: pd.DataFrame, cluster_column: str = "leiden") -> pd.DataFrame:
    summary = (
        obs.groupby([cluster_column, "cluster_cell_type"])
        .size()
        .reset_index(name="n_cells")
        .query("n_cells > 0")
        .sort_values([cluster_column, "n_cells"], ascending=[True, False])
    )
    return summary


def top_marker_table(adata, n_markers: int = 10) -> pd.DataFrame:
    ranked = adata.uns["rank_genes_groups"]
    groups = ranked["names"].dtype.names
    records = []
    for group in groups:
        names = ranked["names"][group][:n_markers]
        scores = ranked["scores"][group][:n_markers]
        pvals_adj = ranked["pvals_adj"][group][:n_markers]
        logfc = ranked["logfoldchanges"][group][:n_markers]
        for gene, score, pval, fc in zip(names, scores, pvals_adj, logfc):
            records.append(
                {
                    "cluster": group,
                    "gene": gene,
                    "score": float(score),
                    "adj_p_value": float(pval),
                    "logfoldchange": float(fc),
                }
            )
    return pd.DataFrame(records)


def write_summary_report(
    destination: Path,
    dataset_name: str,
    cells_before: int,
    cells_after: int,
    genes_before: int,
    genes_after: int,
    cluster_summary: pd.DataFrame,
    marker_table: pd.DataFrame,
) -> None:
    lines = [
        "# Single-cell RNA-seq Analysis Summary",
        "",
        f"- Dataset: `{dataset_name}`",
        f"- Cells before QC: {cells_before}",
        f"- Cells after QC: {cells_after}",
        f"- Genes before QC: {genes_before}",
        f"- Genes after QC: {genes_after}",
        f"- Leiden clusters identified: {cluster_summary['leiden'].nunique()}",
        "",
        "## Cluster Summary",
        "",
    ]

    for _, row in cluster_summary.iterrows():
        lines.append(
            f"- Cluster {row['leiden']}: {row['cluster_cell_type']} ({int(row['n_cells'])} cells)"
        )

    lines.extend(["", "## Representative Marker Genes", ""])
    for cluster, frame in marker_table.groupby("cluster"):
        top_genes = ", ".join(frame["gene"].head(5).tolist())
        lines.append(f"- Cluster {cluster}: {top_genes}")

    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def select_heatmap_markers(marker_table: pd.DataFrame, top_n_per_cluster: int = 3) -> list[str]:
    markers = []
    for _, frame in marker_table.groupby("cluster"):
        markers.extend(frame["gene"].head(top_n_per_cluster).tolist())
    return list(dict.fromkeys(markers))


def percent_mt_label(obs: pd.DataFrame) -> str:
    for candidate in ("pct_counts_mt", "pct_counts_mito"):
        if candidate in obs.columns:
            return candidate
    raise ValueError("Mitochondrial QC column not found in obs.")

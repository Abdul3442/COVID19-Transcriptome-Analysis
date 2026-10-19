#!/usr/bin/env python3
"""
Transcriptome Analysis of COVID-19: Step 6 - Publication-Quality Visualizations
Author: Abdul Rehman
Dataset: NCBI GEO GSE147507

This script generates publication-grade scientific figures:
1. Fig 1: Sequencing depth and detected gene distributions (QC)
2. Fig 2: 2D Principal Component Analysis (PCA) biplot with confidence clustering
3. Fig 3: Sample-to-sample Pearson correlation heatmap
4. Fig 4: Volcano plot of Differential Expression with labeled hallmark COVID-19 biomarkers
5. Fig 5: Hierarchically clustered z-score heatmap of top differentially expressed genes
6. Fig 6: Pathway enrichment bubble / dotplot of top immunological pathways
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

# Plotting style configuration
plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["axes.edgecolor"] = "#2c3e50"

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROCESSED_DATA_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
QC_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "01_qc")
DEA_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "04_dea")
ENRICHMENT_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "05_enrichment")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "results", "06_figures")


def log(msg: str):
    print(f"[06_GENERATE_FIGURES] {msg}", flush=True)


def plot_fig1_qc_library_sizes(cohort_key: str = "clinical_lung"):
    """Figure 1: Barplot of library sizes and detected genes across samples."""
    qc_file = os.path.join(QC_RESULTS_DIR, "qc_metrics_summary.tsv")
    if not os.path.exists(qc_file):
        return

    df_qc = pd.read_csv(qc_file, sep="\t")
    cohort_qc = df_qc[df_qc["cohort"] == cohort_key].copy()
    if cohort_qc.empty:
        cohort_qc = df_qc

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # Color palette
    colors = ["#e74c3c" if "COVID" in str(g) else "#3498db" for g in cohort_qc["group"]]

    # Panel A: Sequencing Depth
    ax1.bar(
        range(len(cohort_qc)),
        cohort_qc["total_reads"] / 1e6,
        color=colors,
        edgecolor="#2c3e50",
        alpha=0.85,
        width=0.6
    )
    ax1.set_title("A. Sequencing Depth (Library Size)", fontsize=13, fontweight="bold", pad=12)
    ax1.set_ylabel("Total Mapped Reads (Millions)", fontsize=11)
    ax1.set_xticks(range(len(cohort_qc)))
    clean_ticks = [s.replace("Series15_", "").replace("Series1_", "").replace("Series7_", "") for s in cohort_qc["sample_id"]]
    ax1.set_xticklabels(clean_ticks, rotation=35, ha="right", fontsize=9)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)

    # Panel B: Genes Detected
    ax2.bar(
        range(len(cohort_qc)),
        cohort_qc["genes_detected_ge10"],
        color=colors,
        edgecolor="#2c3e50",
        alpha=0.85,
        width=0.6
    )
    ax2.set_title("B. Detected Genes (Counts >= 10)", fontsize=13, fontweight="bold", pad=12)
    ax2.set_ylabel("Number of Genes", fontsize=11)
    ax2.set_xticks(range(len(cohort_qc)))
    ax2.set_xticklabels(clean_ticks, rotation=35, ha="right", fontsize=9)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#e74c3c", edgecolor="#2c3e50", label="COVID-19"),
        Patch(facecolor="#3498db", edgecolor="#2c3e50", label="Control")
    ]
    ax1.legend(handles=legend_elements, loc="upper right", frameon=True)

    out_path = os.path.join(FIGURES_DIR, f"{cohort_key}_fig1_qc_library_sizes.png")
    fig.savefig(out_path)
    plt.close(fig)
    log(f"Saved Figure 1 (QC Library Sizes) -> {out_path}")


def plot_fig2_pca(cohort_key: str = "clinical_lung"):
    """Figure 2: 2D Principal Component Analysis Biplot."""
    coords_file = os.path.join(QC_RESULTS_DIR, f"{cohort_key}_pca_coordinates.tsv")
    var_file = os.path.join(QC_RESULTS_DIR, f"{cohort_key}_pca_variance.tsv")
    if not os.path.exists(coords_file) or not os.path.exists(var_file):
        return

    coords_df = pd.read_csv(coords_file, sep="\t")
    var_df = pd.read_csv(var_file, sep="\t")

    pc1_var = var_df.loc[var_df["PC"] == "PC1", "explained_variance_percentage"].values[0]
    pc2_var = var_df.loc[var_df["PC"] == "PC2", "explained_variance_percentage"].values[0]

    fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)

    covid_mask = coords_df["group"] == "COVID19"
    ctrl_mask = coords_df["group"] == "Control"

    ax.scatter(
        coords_df.loc[ctrl_mask, "PC1"],
        coords_df.loc[ctrl_mask, "PC2"],
        c="#3498db",
        s=160,
        edgecolors="#1a5276",
        linewidth=1.8,
        label="Healthy Control",
        zorder=4
    )

    ax.scatter(
        coords_df.loc[covid_mask, "PC1"],
        coords_df.loc[covid_mask, "PC2"],
        c="#e74c3c",
        s=160,
        edgecolors="#78281f",
        linewidth=1.8,
        label="COVID-19 Patient",
        zorder=4
    )

    # Label points
    for _, row in coords_df.iterrows():
        sample_label = row["sample_id"].replace("Series15_", "").replace("Series1_", "").replace("Series7_", "")
        ax.annotate(
            sample_label,
            (row["PC1"], row["PC2"]),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9.5,
            fontweight="semibold",
            color="#2c3e50"
        )

    ax.axhline(0, color="#bdc3c7", linestyle="--", linewidth=1, zorder=1)
    ax.axvline(0, color="#bdc3c7", linestyle="--", linewidth=1, zorder=1)

    ax.set_title(f"Principal Component Analysis ({cohort_key.replace('_', ' ').title()})", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel(f"Principal Component 1 ({pc1_var:.1f}% Variance Explained)", fontsize=11, fontweight="semibold")
    ax.set_ylabel(f"Principal Component 2 ({pc2_var:.1f}% Variance Explained)", fontsize=11, fontweight="semibold")
    ax.legend(loc="best", frameon=True, fancybox=True, shadow=True, fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.6)

    out_path = os.path.join(FIGURES_DIR, f"{cohort_key}_fig2_pca_clustering.png")
    fig.savefig(out_path)
    plt.close(fig)
    log(f"Saved Figure 2 (PCA Clustering) -> {out_path}")


def plot_fig3_sample_correlation(cohort_key: str = "clinical_lung"):
    """Figure 3: Sample-to-sample Pearson correlation matrix heatmap."""
    corr_file = os.path.join(QC_RESULTS_DIR, f"{cohort_key}_sample_correlation.tsv")
    if not os.path.exists(corr_file):
        return

    corr_df = pd.read_csv(corr_file, sep="\t", index_col=0)
    clean_labels = [c.replace("Series15_", "").replace("Series1_", "").replace("Series7_", "") for c in corr_df.columns]
    corr_df.columns = clean_labels
    corr_df.index = clean_labels

    fig, ax = plt.subplots(figsize=(7.5, 6.5), constrained_layout=True)
    sns.heatmap(
        corr_df,
        annot=True,
        fmt=".3f",
        cmap="vlag",
        vmin=0.7,
        vmax=1.0,
        square=True,
        linewidths=1.0,
        cbar_kws={"label": "Pearson Correlation (r)"},
        ax=ax
    )
    ax.set_title(f"Sample-to-Sample Correlation Heatmap\n({cohort_key.replace('_', ' ').title()})", fontsize=13, fontweight="bold", pad=12)

    out_path = os.path.join(FIGURES_DIR, f"{cohort_key}_fig3_sample_correlation_heatmap.png")
    fig.savefig(out_path)
    plt.close(fig)
    log(f"Saved Figure 3 (Sample Correlation Heatmap) -> {out_path}")


def plot_fig4_volcano(cohort_key: str = "clinical_lung"):
    """Figure 4: Volcano plot of Differential Expression with key COVID-19 biomarkers labeled."""
    dea_file = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_dea_all_genes.tsv")
    if not os.path.exists(dea_file):
        return

    dea_df = pd.read_csv(dea_file, sep="\t", index_col=0)

    # Compute -log10(padj)
    clamped_padj = np.clip(dea_df["padj"].values, 1e-30, 1.0)
    neg_log10_padj = -np.log10(clamped_padj)
    dea_df["neg_log10_padj"] = neg_log10_padj

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)

    # Non-significant
    non_sig = dea_df[dea_df["status"] == "Not_Significant"]
    ax.scatter(non_sig["log2_fold_change"], non_sig["neg_log10_padj"], c="#bdc3c7", s=18, alpha=0.5, label="Not Significant", zorder=2)

    # Upregulated
    up = dea_df[dea_df["status"] == "Significant_Upregulated"]
    ax.scatter(up["log2_fold_change"], up["neg_log10_padj"], c="#e74c3c", s=28, alpha=0.8, label=f"Upregulated ({len(up)})", zorder=3)

    # Downregulated
    down = dea_df[dea_df["status"] == "Significant_Downregulated"]
    ax.scatter(down["log2_fold_change"], down["neg_log10_padj"], c="#3498db", s=28, alpha=0.8, label=f"Downregulated ({len(down)})", zorder=3)

    # Cutoff threshold lines
    ax.axvline(1.0, color="#7f8c8d", linestyle="--", linewidth=1.2, alpha=0.7)
    ax.axvline(-1.0, color="#7f8c8d", linestyle="--", linewidth=1.2, alpha=0.7)
    ax.axhline(-np.log10(0.05), color="#7f8c8d", linestyle="--", linewidth=1.2, alpha=0.7)

    # Specific hallmark genes to label
    key_markers = [
        "CXCL8", "IL6", "CXCL1", "CXCL2", "CCL2", "CCL20", "TNF", "NFKBIA",
        "IFI6", "ISG15", "OAS1", "OAS2", "IFIT1", "IFIT2", "MX1", "STAT1", "ACE2", "S100A8", "S100A9"
    ]
    labeled_genes = set()
    for gene in key_markers:
        if gene in dea_df.index and dea_df.loc[gene, "status"] != "Not_Significant":
            row = dea_df.loc[gene]
            ax.annotate(
                gene,
                (row["log2_fold_change"], row["neg_log10_padj"]),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=9,
                fontweight="bold",
                color="#800020",
                arrowprops=dict(arrowstyle="->", color="#800020", lw=0.8)
            )
            labeled_genes.add(gene)

    # Also annotate top 5 most extreme upregulated genes if not already labeled
    top_up_candidates = up.sort_values(by="neg_log10_padj", ascending=False).head(8).index
    for gene in top_up_candidates:
        if gene not in labeled_genes:
            row = dea_df.loc[gene]
            ax.annotate(
                gene,
                (row["log2_fold_change"], row["neg_log10_padj"]),
                xytext=(8, 8),
                textcoords="offset points",
                fontsize=8.5,
                fontweight="semibold",
                color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.6)
            )
            labeled_genes.add(gene)

    ax.set_title(f"Differential Expression Volcano Plot\nCOVID-19 vs Control ({cohort_key.replace('_', ' ').title()})", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel(r"$\log_2$ Fold Change (COVID-19 / Control)", fontsize=11, fontweight="semibold")
    ax.set_ylabel(r"$-\log_{10}$ (Adjusted $p$-value / FDR)", fontsize=11, fontweight="semibold")
    ax.legend(loc="upper left", frameon=True, fancybox=True, shadow=True, fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.5)

    out_path = os.path.join(FIGURES_DIR, f"{cohort_key}_fig4_volcano_plot_degs.png")
    fig.savefig(out_path)
    plt.close(fig)
    log(f"Saved Figure 4 (Volcano Plot) -> {out_path}")


def plot_fig5_heatmap_top_degs(cohort_key: str = "clinical_lung", top_n: int = 35):
    """Figure 5: Hierarchically clustered z-score heatmap of top differentially expressed genes."""
    log2_file = os.path.join(PROCESSED_DATA_DIR, f"{cohort_key}_log2_normalized_counts.tsv")
    dea_file = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_dea_all_genes.tsv")
    meta_file = os.path.join(PROCESSED_DATA_DIR, f"{cohort_key}_metadata.tsv")

    if not os.path.exists(log2_file) or not os.path.exists(dea_file) or not os.path.exists(meta_file):
        return

    log2_df = pd.read_csv(log2_file, sep="\t", index_col=0)
    dea_df = pd.read_csv(dea_file, sep="\t", index_col=0)
    meta_df = pd.read_csv(meta_file, sep="\t").set_index("sample_id")

    # Select top DEGs by adjusted p-value
    sig_degs = dea_df[dea_df["status"] != "Not_Significant"].sort_values(by="padj")
    if sig_degs.empty:
        # If no significant DEGs meet the threshold, take top by absolute log2FC
        sig_degs = dea_df.sort_values(by="log2_fold_change", key=abs, ascending=False)

    selected_genes = sig_degs.head(min(top_n, len(sig_degs))).index.tolist()
    if len(selected_genes) < 2:
        return

    sub_matrix = log2_df.loc[selected_genes].copy()

    # Z-score normalization per row (gene): (x - mean) / std
    z_matrix = sub_matrix.apply(lambda row: (row - row.mean()) / (row.std() if row.std() > 0 else 1.0), axis=1)

    # Clean sample names
    z_matrix.columns = [c.replace("Series15_", "").replace("Series1_", "").replace("Series7_", "") for c in z_matrix.columns]

    # Create clustermap
    group_map = meta_df["group"].to_dict()
    clean_sample_groups = [group_map.get(s, "Unknown") for s in meta_df.index]
    col_colors = ["#e74c3c" if "COVID" in str(g) else "#3498db" for g in clean_sample_groups]

    try:
        g = sns.clustermap(
            z_matrix,
            cmap="coolwarm",
            center=0,
            vmin=-2.0,
            vmax=2.0,
            col_colors=col_colors if len(col_colors) == z_matrix.shape[1] else None,
            col_cluster=True,
            row_cluster=True,
            figsize=(8.5, 10.5),
            cbar_kws={"label": "Expression Z-Score"},
            dendrogram_ratio=(0.15, 0.15),
            linewidths=0.5,
            linecolor="#f0f0f0"
        )
        g.ax_heatmap.set_title(f"Top Differentially Expressed Genes\n({cohort_key.replace('_', ' ').title()})", fontsize=13, fontweight="bold", pad=20)
        g.ax_heatmap.set_xlabel("Samples", fontsize=10, fontweight="semibold")
        g.ax_heatmap.set_ylabel("Genes", fontsize=10, fontweight="semibold")

        out_path = os.path.join(FIGURES_DIR, f"{cohort_key}_fig5_top_degs_heatmap.png")
        g.savefig(out_path)
        plt.close("all")
        log(f"Saved Figure 5 (Clustered Heatmap) -> {out_path}")
    except Exception as e:
        log(f"Heatmap plotting exception for {cohort_key}: {e}")


def plot_fig6_pathway_bubble(cohort_key: str = "clinical_lung"):
    """Figure 6: Pathway enrichment bubble / dotplot."""
    up_ora_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_upregulated_pathways.tsv")
    all_ora_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_all_degs_pathways.tsv")
    enrichr_go_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_enrichr_go_bp.tsv")

    enrich_df = None
    if os.path.exists(up_ora_file) and os.path.getsize(up_ora_file) > 10:
        enrich_df = pd.read_csv(up_ora_file, sep="\t")
    elif os.path.exists(all_ora_file) and os.path.getsize(all_ora_file) > 10:
        enrich_df = pd.read_csv(all_ora_file, sep="\t")
    elif os.path.exists(enrichr_go_file) and os.path.getsize(enrichr_go_file) > 10:
        enrich_df = pd.read_csv(enrichr_go_file, sep="\t")

    if enrich_df is None or enrich_df.empty:
        return

    # Filter to top pathways
    top_pathways = enrich_df.head(10).copy()
    top_pathways = top_pathways.sort_values(by="padj", ascending=False)

    clamped_padj = np.clip(top_pathways["padj"].values, 1e-15, 1.0)
    top_pathways["neg_log10_padj"] = -np.log10(clamped_padj)

    fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)

    sizes = top_pathways["overlap_count"] * 45 + 50
    scatter = ax.scatter(
        top_pathways["neg_log10_padj"],
        range(len(top_pathways)),
        s=sizes,
        c=top_pathways["neg_log10_padj"],
        cmap="YlOrRd",
        edgecolors="#2c3e50",
        linewidth=1.2,
        alpha=0.9
    )

    ax.set_yticks(range(len(top_pathways)))
    clean_names = [p[:50] + "..." if len(p) > 50 else p for p in top_pathways["pathway"]]
    ax.set_yticklabels(clean_names, fontsize=10, fontweight="semibold")

    ax.set_xlabel(r"$-\log_{10}$ (Adjusted $p$-value)", fontsize=11, fontweight="semibold")
    ax.set_title(f"Top Enriched Biological Pathways in COVID-19\n({cohort_key.replace('_', ' ').title()})", fontsize=13, fontweight="bold", pad=12)
    ax.grid(True, linestyle=":", alpha=0.6)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label(r"$-\log_{10}$ (FDR)", fontsize=10, fontweight="semibold")

    from matplotlib.lines import Line2D
    legend_sizes = [5, 10, 15]
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', label=f"{s} genes",
               markerfacecolor='#e67e22', markeredgecolor='#2c3e50', markersize=np.sqrt(s*45+50)/2)
        for s in legend_sizes if s <= top_pathways["overlap_count"].max() + 5
    ]
    if legend_handles:
        ax.legend(handles=legend_handles, title="Gene Count", loc="lower right", frameon=True)

    out_path = os.path.join(FIGURES_DIR, f"{cohort_key}_fig6_pathway_enrichment_dotplot.png")
    fig.savefig(out_path)
    plt.close(fig)
    log(f"Saved Figure 6 (Pathway Enrichment Dotplot) -> {out_path}")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    cohorts = ["clinical_lung", "primary_nhbe", "calu3_cells"]
    for c in cohorts:
        log(f"Generating figures for cohort '{c}'...")
        plot_fig1_qc_library_sizes(c)
        plot_fig2_pca(c)
        plot_fig3_sample_correlation(c)
        plot_fig4_volcano(c)
        plot_fig5_heatmap_top_degs(c)
        plot_fig6_pathway_bubble(c)
    log("Step 6 (Publication-Quality Visualizations) completed successfully!")


if __name__ == "__main__":
    main()

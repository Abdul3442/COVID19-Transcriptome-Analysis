#!/usr/bin/env python3
"""
Transcriptome Analysis of COVID-19: Step 4 - Differential Expression Analysis (DEA)
Author: Abdul Rehman
Dataset: NCBI GEO GSE147507

This script:
1. Performs rigorous differential expression testing between COVID-19 and Control.
2. Estimates base mean, condition-specific means, log2 fold change (log2FC), and standard error.
3. Calculates statistical significance via empirical-moderated t / Wald statistics.
4. Corrects for multiple testing using the Benjamini-Hochberg False Discovery Rate (FDR).
5. Categorizes genes (Upregulated, Downregulated, Not Significant).
6. Identifies and isolates key SARS-CoV-2 signature genes (cytokines, ISGs, NF-kB targets).
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROCESSED_DATA_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
DEA_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "04_dea")


def log(msg: str):
    print(f"[04_DIFFERENTIAL_EXPRESSION] {msg}", flush=True)


def benjamini_hochberg_correction(p_values: np.ndarray) -> np.ndarray:
    """
    Computes Benjamini-Hochberg (BH) adjusted p-values (False Discovery Rate, FDR).
    Maintains monotonicity so that higher nominal p-values cannot have lower adjusted p-values.
    """
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.array([])

    order = np.argsort(p_values)
    ranked_p = p_values[order]
    
    # Calculate step-up multipliers: (p * n) / rank
    ranks = np.arange(1, n + 1)
    adj_p = ranked_p * (n / ranks)
    
    # Enforce monotonicity backwards: adj_p[i] = min(adj_p[i], adj_p[i+1])
    adj_p = np.minimum.accumulate(adj_p[::-1])[::-1]
    # Bound by 1.0
    adj_p = np.clip(adj_p, 0.0, 1.0)
    
    # Re-order to original index
    original_adj_p = np.empty(n, dtype=float)
    original_adj_p[order] = adj_p
    return original_adj_p


def run_differential_expression(norm_counts_df: pd.DataFrame, 
                                log2_counts_df: pd.DataFrame, 
                                meta_df: pd.DataFrame,
                                cohort_name: str,
                                lfc_threshold: float = 1.0,
                                fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Performs differential gene expression testing between COVID-19 and Control.
    """
    covid_samples = meta_df[meta_df["group"] == "COVID19"]["sample_id"].tolist()
    control_samples = meta_df[meta_df["group"] == "Control"]["sample_id"].tolist()

    n_covid = len(covid_samples)
    n_ctrl = len(control_samples)
    log(f"Cohort '{cohort_name}': Comparing {n_covid} COVID-19 vs {n_ctrl} Control samples...")

    results = []
    genes = norm_counts_df.index.tolist()

    norm_ctrl = norm_counts_df[control_samples].values
    norm_covid = norm_counts_df[covid_samples].values
    
    log2_ctrl = log2_counts_df[control_samples].values
    log2_covid = log2_counts_df[covid_samples].values

    # Overall base mean
    base_means = np.mean(norm_counts_df.values, axis=1)
    mean_ctrl = np.mean(norm_ctrl, axis=1)
    mean_covid = np.mean(norm_covid, axis=1)

    # Log2 fold change computed from log2(norm + 1) means
    mean_log2_ctrl = np.mean(log2_ctrl, axis=1)
    mean_log2_covid = np.mean(log2_covid, axis=1)
    log2_fc = mean_log2_covid - mean_log2_ctrl

    # Variance and Standard Error estimation across replicates
    var_ctrl = np.var(log2_ctrl, axis=1, ddof=1) if n_ctrl > 1 else np.zeros(len(genes))
    var_covid = np.var(log2_covid, axis=1, ddof=1) if n_covid > 1 else np.zeros(len(genes))
    
    # Degrees of freedom
    df = (n_ctrl - 1) + (n_covid - 1)
    if df < 1:
        df = 1

    # Pooled variance
    pooled_var = ((n_ctrl - 1) * var_ctrl + (n_covid - 1) * var_covid) / max(1, df)
    
    # Variance stabilization floor (prevents artificially zero variance in identical small counts)
    min_variance_floor = np.percentile(pooled_var[pooled_var > 0], 5) if np.any(pooled_var > 0) else 1e-4
    pooled_var_moderated = np.maximum(pooled_var, min_variance_floor)

    # Standard error of the difference
    se = np.sqrt(pooled_var_moderated * (1.0 / n_ctrl + 1.0 / n_covid))

    # Test statistic (t-statistic)
    t_stat = log2_fc / se

    # Two-tailed p-values from t-distribution
    p_values = 2.0 * (1.0 - stats.t.cdf(np.abs(t_stat), df=df))
    # Safeguard against NaN or out of bounds
    p_values = np.nan_to_num(p_values, nan=1.0, posinf=1.0, neginf=1.0)
    p_values = np.clip(p_values, 1e-300, 1.0)

    # False Discovery Rate adjustment
    fdr_values = benjamini_hochberg_correction(p_values)

    # Compile into DataFrame
    dea_df = pd.DataFrame({
        "gene_symbol": genes,
        "base_mean": np.round(base_means, 2),
        "mean_control": np.round(mean_ctrl, 2),
        "mean_covid19": np.round(mean_covid, 2),
        "log2_fold_change": np.round(log2_fc, 4),
        "std_error": np.round(se, 4),
        "t_statistic": np.round(t_stat, 4),
        "pvalue": p_values,
        "padj": fdr_values
    }).set_index("gene_symbol")

    # Classify regulation status
    conditions = [
        (dea_df["padj"] < fdr_threshold) & (dea_df["log2_fold_change"] >= lfc_threshold),
        (dea_df["padj"] < fdr_threshold) & (dea_df["log2_fold_change"] <= -lfc_threshold),
    ]
    choices = ["Significant_Upregulated", "Significant_Downregulated"]
    dea_df["status"] = np.select(conditions, choices, default="Not_Significant")

    return dea_df


def process_cohort_dea(cohort_key: str):
    norm_file = os.path.join(PROCESSED_DATA_DIR, f"{cohort_key}_normalized_counts.tsv")
    log2_file = os.path.join(PROCESSED_DATA_DIR, f"{cohort_key}_log2_normalized_counts.tsv")
    meta_file = os.path.join(PROCESSED_DATA_DIR, f"{cohort_key}_metadata.tsv")

    if not os.path.exists(norm_file) or not os.path.exists(meta_file):
        log(f"Normalized files not found for {cohort_key}. Skipping DEA.")
        return

    norm_counts = pd.read_csv(norm_file, sep="\t", index_col=0)
    log2_counts = pd.read_csv(log2_file, sep="\t", index_col=0)
    meta_df = pd.read_csv(meta_file, sep="\t")

    dea_df = run_differential_expression(
        norm_counts_df=norm_counts,
        log2_counts_df=log2_counts,
        meta_df=meta_df,
        cohort_name=cohort_key,
        lfc_threshold=1.0,
        fdr_threshold=0.05
    )

    # Summary statistics
    n_total = len(dea_df)
    n_up = (dea_df["status"] == "Significant_Upregulated").sum()
    n_down = (dea_df["status"] == "Significant_Downregulated").sum()
    n_not_sig = (dea_df["status"] == "Not_Significant").sum()

    log(f"Cohort '{cohort_key}' Results: Total={n_total}, Upregulated={n_up} ({n_up/n_total*100:.1f}%), "
        f"Downregulated={n_down} ({n_down/n_total*100:.1f}%), NotSig={n_not_sig}")

    # Top tables
    up_genes = dea_df[dea_df["status"] == "Significant_Upregulated"].sort_values(
        by=["padj", "log2_fold_change"], ascending=[True, False]
    )
    down_genes = dea_df[dea_df["status"] == "Significant_Downregulated"].sort_values(
        by=["padj", "log2_fold_change"], ascending=[True, True]
    )

    # Save outputs
    out_all = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_dea_all_genes.tsv")
    out_up = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_top_upregulated_genes.tsv")
    out_down = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_top_downregulated_genes.tsv")
    out_summary = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_dea_summary.tsv")

    dea_df.to_csv(out_all, sep="\t")
    up_genes.head(50).to_csv(out_up, sep="\t")
    down_genes.head(50).to_csv(out_down, sep="\t")

    summary_df = pd.DataFrame([{
        "cohort": cohort_key,
        "total_genes_tested": n_total,
        "significant_upregulated": n_up,
        "significant_downregulated": n_down,
        "not_significant": n_not_sig,
        "lfc_cutoff": 1.0,
        "fdr_cutoff": 0.05
    }])
    summary_df.to_csv(out_summary, sep="\t", index=False)

    # Highlight canonical COVID-19 hallmark genes
    hallmark_genes = [
        "IL6", "CXCL8", "CXCL1", "CXCL2", "CCL2", "CCL20", "TNF", "NFKBIA",
        "IFIT1", "IFIT2", "IFIT3", "ISG15", "OAS1", "OAS2", "IFI6", "MX1", "STAT1", "ACE2"
    ]
    present_hallmarks = [g for g in hallmark_genes if g in dea_df.index]
    if present_hallmarks:
        hallmark_df = dea_df.loc[present_hallmarks][["base_mean", "log2_fold_change", "pvalue", "padj", "status"]]
        hallmark_out = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_hallmark_covid_genes.tsv")
        hallmark_df.to_csv(hallmark_out, sep="\t")
        log(f"Found {len(present_hallmarks)} COVID-19 hallmark genes. Saved to: {hallmark_out}")


def main():
    os.makedirs(DEA_RESULTS_DIR, exist_ok=True)
    cohorts = ["clinical_lung", "primary_nhbe", "calu3_cells"]
    for c in cohorts:
        process_cohort_dea(c)
    log("Step 4 (Differential Expression Analysis) completed successfully!")


if __name__ == "__main__":
    main()

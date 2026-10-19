#!/usr/bin/env python3
"""
Transcriptome Analysis of COVID-19: Step 5 - Functional Pathway Enrichment Analysis
Author: Abdul Rehman
Dataset: NCBI GEO GSE147507

This script:
1. Performs Over-Representation Analysis (ORA) on significantly upregulated and downregulated genes.
2. Evaluates Gene Ontology (GO Biological Process), KEGG pathways, and Reactome pathways.
3. Implements dual-mode enrichment:
   - Remote: Public Enrichr REST API (maayanlab.cloud) for comprehensive library scanning.
   - Local: Fisher's exact / hypergeometric test with Benjamini-Hochberg FDR correction on curated immune gene sets.
4. Exports structured enrichment tables with overlap counts, odds ratios, p-values, and FDR adjustments.
"""

import os
import sys
import json
import time
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEA_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "04_dea")
ENRICHMENT_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "05_enrichment")


def log(msg: str):
    print(f"[05_PATHWAY_ENRICHMENT] {msg}", flush=True)


# Curated Hallmark COVID-19 & Innate Immune Gene Sets for guaranteed local testing
CURATED_PATHWAYS = {
    "GO: Defense response to virus": [
        "OAS1", "OAS2", "OAS3", "OASL", "MX1", "MX2", "ISG15", "IFIT1", "IFIT2", "IFIT3", 
        "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L", "STAT1", "STAT2", "IRF7", "IRF9", "RSAD2", 
        "DDX58", "IFIH1", "EIF2AK2", "BST2", "HERC5", "XAF1", "USP18", "TRIM22"
    ],
    "GO: Cytokine-mediated signaling pathway": [
        "IL6", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CCL2", "CCL3", "CCL4", "CCL20", "TNF",
        "IL1B", "IL1A", "IL1RN", "IL10", "SOCS3", "NFKBIA", "NFKBIZ", "JUN", "FOS", "STAT3",
        "TNFAIP3", "PTGS2", "ICAM1", "CSF2", "CSF3", "IL18"
    ],
    "GO: Response to type I interferon": [
        "IFIT1", "IFIT2", "IFIT3", "ISG15", "OAS1", "OAS2", "OAS3", "MX1", "MX2", "IFI6",
        "IFI35", "STAT1", "STAT2", "IRF9", "USP18", "RSAD2", "BST2", "HERC5", "XAF1", "SAMD9"
    ],
    "GO: Regulation of inflammatory response": [
        "IL6", "CXCL8", "TNF", "IL1B", "NFKBIA", "NFKB1", "RELA", "TLR2", "TLR4", "MYD88",
        "NOD2", "NLRP3", "CASP1", "IL18", "PTGS2", "SERPINE1", "CXCL10", "CCL2", "CCL5"
    ],
    "GO: Neutrophil degranulation / activation": [
        "S100A8", "S100A9", "S100A12", "LCN2", "MMP9", "CAMP", "DEFA1", "DEFA3", "LTF",
        "ELANE", "MPO", "BPI", "CTSG", "AZU1", "PRTN3", "CD14", "FCGR3B", "ALPL"
    ],
    "KEGG: Coronavirus disease - COVID-19": [
        "ACE2", "TMPRSS2", "STAT1", "STAT3", "JAK1", "JAK2", "IL6", "CXCL8", "TNF", "IL1B",
        "NFKBIA", "NFKB1", "RELA", "OAS1", "OAS2", "MX1", "ISG15", "IFIT1", "JUN", "FOS",
        "CXCL10", "CCL2", "CCL3", "CCL5", "CASP1", "NLRP3", "IL18"
    ],
    "KEGG: NF-kappa B signaling pathway": [
        "NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA", "NFKBIB", "NFKBIZ", "TNFAIP3", "TRAF2",
        "TRAF6", "CHUK", "IKBKB", "IKBKG", "RIPK1", "TAB1", "TAB2", "MAP3K7", "CXCL8", "IL1B", "TNF"
    ],
    "KEGG: TNF signaling pathway": [
        "TNF", "TNFRSF1A", "TNFRSF1B", "TRADD", "FADD", "TRAF2", "CASP8", "CASP3", "RIPK1",
        "NFKBIA", "NFKB1", "RELA", "MAPK8", "JUN", "FOS", "CXCL1", "CXCL2", "CXCL8", "IL6", "CCL2", "PTGS2"
    ],
    "KEGG: Cytokine-cytokine receptor interaction": [
        "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CXCL10", "CCL2", "CCL3", "CCL4", "CCL5", "CCL20",
        "IL6", "IL6R", "IL1B", "IL1R1", "TNF", "TNFRSF1A", "IL10", "IL10RA", "CSF2", "CSF3", "IFNAR1"
    ],
    "Reactome: Interleukin-6 signaling": [
        "IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "TYK2", "STAT3", "SOCS3", "PTPN11", "MAPK1", "MAPK3"
    ],
    "Reactome: Chemokine receptors bind chemokines": [
        "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CXCR1", "CXCR2", "CCL2", "CCR2", "CCL20", "CCR6", "GNAI1"
    ],
}


def benjamini_hochberg(p_vals: np.ndarray) -> np.ndarray:
    n = len(p_vals)
    if n == 0:
        return np.array([])
    order = np.argsort(p_vals)
    ranked = p_vals[order]
    adj = ranked * (n / np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    res = np.empty(n)
    res[order] = adj
    return res


def run_local_ora(deg_genes: list, universe_genes: list) -> pd.DataFrame:
    """
    Performs Fisher's Exact test (Hypergeometric distribution) for 2x2 contingency table:
    [ [genes_in_set_and_deg,     genes_in_set_not_deg],
      [genes_not_set_and_deg, genes_not_set_not_deg] ]
    """
    deg_set = set(deg_genes)
    universe_set = set(universe_genes)
    N = len(universe_set)
    n = len(deg_set)

    records = []
    for term, term_genes in CURATED_PATHWAYS.items():
        K_genes = [g for g in term_genes if g in universe_set]
        K = len(K_genes)
        if K < 3:
            continue
        overlap_genes = [g for g in term_genes if g in deg_set]
        k = len(overlap_genes)

        # 2x2 contingency table
        table = [
            [k, K - k],
            [n - k, N - (K + n - k)]
        ]
        odds_ratio, p_val = stats.fisher_exact(table, alternative="greater")

        records.append({
            "pathway": term,
            "overlap_count": k,
            "pathway_gene_count": K,
            "deg_query_count": n,
            "odds_ratio": round(odds_ratio, 3) if not np.isinf(odds_ratio) else 999.0,
            "pvalue": p_val,
            "genes": ";".join(sorted(overlap_genes))
        })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    df["padj"] = benjamini_hochberg(df["pvalue"].values)
    df = df.sort_values(by=["padj", "pvalue"], ascending=[True, True])
    return df


def query_enrichr_api(gene_list: list, library_name: str = "GO_Biological_Process_2023") -> pd.DataFrame:
    """Queries official Enrichr REST API for comprehensive library testing."""
    if len(gene_list) < 5:
        return pd.DataFrame()

    post_url = "https://maayanlab.cloud/Enrichr/addList"
    genes_str = "\n".join(gene_list)
    payload = urllib.parse.urlencode({
        "list": genes_str,
        "description": "COVID19_DEA_Gene_List"
    }).encode("utf-8")

    try:
        req = urllib.request.Request(post_url, data=payload, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        user_list_id = data.get("userListId")
        if not user_list_id:
            return pd.DataFrame()

        get_url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={library_name}"
        time.sleep(0.5)
        req_get = urllib.request.Request(get_url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req_get, timeout=20) as resp_get:
            res_data = json.loads(resp_get.read().decode("utf-8"))

        entries = res_data.get(library_name, [])
        records = []
        for e in entries:
            # Format: [rank, term, pvalue, zscore, combined_score, genes, adj_pvalue, ...]
            records.append({
                "pathway": e[1],
                "pvalue": float(e[2]),
                "combined_score": float(e[4]),
                "padj": float(e[6]),
                "genes": ";".join(e[5]) if isinstance(e[5], list) else str(e[5]),
                "overlap_count": len(e[5]) if isinstance(e[5], list) else 0
            })
        return pd.DataFrame(records)
    except Exception as exc:
        log(f"Enrichr API remote query encountered: {exc}. Proceeding with local ORA.")
        return pd.DataFrame()


def process_cohort_enrichment(cohort_key: str):
    dea_file = os.path.join(DEA_RESULTS_DIR, f"{cohort_key}_dea_all_genes.tsv")
    if not os.path.exists(dea_file):
        log(f"DEA file not found for {cohort_key}. Skipping enrichment.")
        return

    dea_df = pd.read_csv(dea_file, sep="\t", index_col=0)
    universe_genes = dea_df.index.tolist()

    up_genes = dea_df[dea_df["status"] == "Significant_Upregulated"].index.tolist()
    down_genes = dea_df[dea_df["status"] == "Significant_Downregulated"].index.tolist()
    all_degs = up_genes + down_genes

    log(f"Cohort '{cohort_key}': Running enrichment on {len(up_genes)} UP genes and {len(down_genes)} DOWN genes...")

    # 1. Local ORA for Upregulated genes
    up_ora = run_local_ora(up_genes, universe_genes)
    if not up_ora.empty:
        up_ora_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_upregulated_pathways.tsv")
        up_ora.to_csv(up_ora_file, sep="\t", index=False)
        log(f"Saved local ORA for UP genes -> {up_ora_file}")

    # 2. Local ORA for All DEGs
    all_ora = run_local_ora(all_degs, universe_genes)
    if not all_ora.empty:
        all_ora_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_all_degs_pathways.tsv")
        all_ora.to_csv(all_ora_file, sep="\t", index=False)

    # 3. Try Remote Enrichr API (GO Biological Process and KEGG)
    enrichr_go = query_enrichr_api(up_genes, library_name="GO_Biological_Process_2023")
    if not enrichr_go.empty:
        go_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_enrichr_go_bp.tsv")
        enrichr_go.head(50).to_csv(go_file, sep="\t", index=False)
        log(f"Saved Enrichr GO Biological Process -> {go_file}")

    enrichr_kegg = query_enrichr_api(up_genes, library_name="KEGG_2021_Human")
    if not enrichr_kegg.empty:
        kegg_file = os.path.join(ENRICHMENT_RESULTS_DIR, f"{cohort_key}_enrichr_kegg.tsv")
        enrichr_kegg.head(50).to_csv(kegg_file, sep="\t", index=False)
        log(f"Saved Enrichr KEGG Pathways -> {kegg_file}")


def main():
    os.makedirs(ENRICHMENT_RESULTS_DIR, exist_ok=True)
    cohorts = ["clinical_lung", "primary_nhbe", "calu3_cells"]
    for c in cohorts:
        process_cohort_enrichment(c)
    log("Step 5 (Pathway Enrichment Analysis) completed successfully!")


if __name__ == "__main__":
    main()

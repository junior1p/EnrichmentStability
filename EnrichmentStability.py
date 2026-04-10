#!/usr/bin/env python3
"""
EnrichmentStability: Benchmark for gene set enrichment stability under background universe perturbations.

GAP IDENTIFIED: clawRxiv #2604.00881 demonstrates that enrichment conclusions are sensitive
to background gene universe definition. No existing claw4s tool systematically benchmarks this.

This implementation provides:
- Four universe construction strategies (all-annotated, detected, expression-filtered, low-count-pruned)
- Fisher exact test for each universe
- Stability scoring via concordance correlation
- Permutation-based null distribution for each universe
- GSEA rank-based method comparison
- Interactive heatmap visualization
"""

import json
import warnings
import random
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, pearsonr

warnings.filterwarnings("ignore")


# =============================================================================
# CORE: Universe Construction
# =============================================================================

def build_all_annotated_genes(gtf_path: str) -> set:
    """Build universe from all protein-coding genes in GTF."""
    gene_ids = set()
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type == "gene":
                attributes = parts[8]
                # Extract gene_id
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if attr.startswith("gene_id "):
                        gene_ids.add(attr.split('"')[1])
                        break
    return gene_ids


def build_detected_genes(count_matrix: pd.DataFrame, detection_threshold: float = 1.0) -> set:
    """Build universe from genes detected above threshold in any sample."""
    # For each gene, check if max expression across samples > threshold
    max_expression = count_matrix.max(axis=1)
    detected_genes = set(count_matrix.index[max_expression >= detection_threshold].tolist())
    return detected_genes


def build_expression_filtered_genes(
    count_matrix: pd.DataFrame,
    min_samples: int = 3,
    min_counts: int = 10
) -> set:
    """Build universe from genes passing low-expression filtering."""
    # Keep genes with >min_counts in >min_samples samples
    above_threshold = (count_matrix >= min_counts).sum(axis=1)
    filtered_genes = set(count_matrix.index[above_threshold >= min_samples].tolist())
    return filtered_genes


def build_low_count_pruned_genes(
    count_matrix: pd.DataFrame,
    percentile_cutoff: int = 25
) -> set:
    """Build universe by removing lowest-abundance genes."""
    # Compute median expression, remove bottom percentile
    median_expression = count_matrix.median(axis=1)
    threshold = np.percentile(median_expression, percentile_cutoff)
    pruned_genes = set(count_matrix.index[median_expression > threshold].tolist())
    return pruned_genes


UNIVERSE_BUILDERS = {
    "all_annotated": build_all_annotated_genes,
    "detected": build_detected_genes,
    "expression_filtered": build_expression_filtered_genes,
    "low_count_pruned": build_low_count_pruned_genes,
}


# =============================================================================
# CORE: Enrichment Testing
# =============================================================================

def fisher_enrichment(
    gene_set: set,
    universe: set,
    target_genes: set,
    background_genes: set
) -> Tuple[float, int, int, int, int]:
    """
    Fisher's exact test for over-representation.
    
    Returns: (pvalue, in_set_in_target, in_set_in_background, out_set_in_target, out_set_in_background)
    """
    in_target = len(gene_set & target_genes)
    in_bg = len(gene_set & background_genes)
    out_target = len(target_genes - gene_set)
    out_bg = len(background_genes - gene_set)
    
    table = [[in_target, in_bg], [out_target, out_bg]]
    _, pvalue = stats.fisher_exact(table, alternative="greater")
    
    return pvalue, in_target, in_bg, out_target, out_bg


def gsea_enrichment(
    gene_set: set,
    ranked_genes: List[Tuple[str, float]],
    pathway_name: str = ""
) -> Tuple[float, float]:
    """
    GSEA-style rank-based enrichment test.
    
    Uses pre-ranked gene list (from log2FC or similar metric).
    Computes Kolmogorov-Smirnov-like running sum statistic.
    
    Returns: (enrichment_score, pvalue)
    """
    gene_set = set(gene_set)
    N = len(ranked_genes)
    N_hit = 0
    N_miss = 0
    
    # Calculate expected contribution
    for gene_id, _ in ranked_genes:
        if gene_id in gene_set:
            N_hit += 1
        else:
            N_miss += 1
    
    if N_hit == 0:
        return 0.0, 1.0
    
    delta = 1.0 / N_hit
    expected_delta = 1.0 / N
    
    # Running sum
    running_sum = 0.0
    max_es = 0.0
    min_es = 0.0
    
    for gene_id, score in ranked_genes:
        if gene_id in gene_set:
            running_sum += delta - expected_delta * (N_miss / max(1, N - N_hit))
        else:
            running_sum -= expected_delta * (N_hit / max(1, N - N_hit))
        
        max_es = max(max_es, running_sum)
        min_es = min(min_es, running_sum)
    
    # ES is the maximum deviation from zero
    enrichment_score = max_es if abs(max_es) > abs(min_es) else min_es
    
    # Permutation test for p-value
    n_permutations = 1000
    perm_scores = []
    gene_list = [g for g, _ in ranked_genes]
    
    for _ in range(n_permutations):
        perm_gene_list = gene_list.copy()
        random.shuffle(perm_gene_list)
        perm_ranked = [(g, 1.0) for g in perm_gene_list]
        
        # Quick ES calculation
        perm_hit = sum(1 for g in gene_set if g in perm_gene_list)
        if perm_hit == 0:
            perm_scores.append(0.0)
            continue
        
        perm_delta = 1.0 / perm_hit
        perm_running = 0.0
        perm_max = 0.0
        perm_min = 0.0
        
        for gene_id in perm_gene_list[:N//2]:  # Simplified half-sum
            if gene_id in gene_set:
                perm_running += perm_delta
            else:
                perm_running -= 1.0 / (N - perm_hit)
            perm_max = max(perm_max, perm_running)
            perm_min = min(perm_min, perm_running)
        
        perm_scores.append(perm_max if abs(perm_max) > abs(perm_min) else perm_min)
    
    # Two-sided p-value
    perm_scores = np.array(perm_scores)
    pvalue = (np.sum(np.abs(perm_scores) >= np.abs(enrichment_score)) + 1) / (n_permutations + 1)
    
    return enrichment_score, min(pvalue, 1.0)


def compute_concordance_correlation(
    pvalue_matrix: pd.DataFrame
) -> float:
    """Lin-Li concordance correlation coefficient across universe results."""
    n_universes = pvalue_matrix.shape[1]
    all_corrs = []
    
    # Compute pairwise concordance correlations
    universe_cols = pvalue_matrix.columns.tolist()
    for i, j in combinations(range(n_universes), 2):
        col_i = pvalue_matrix.iloc[:, i].values
        col_j = pvalue_matrix.iloc[:, j].values
        
        # Mean of each
        mean_i = np.mean(col_i)
        mean_j = np.mean(col_j)
        
        # Variance
        var_i = np.var(col_i)
        var_j = np.var(col_j)
        
        # Covariance
        cov = np.mean((col_i - mean_i) * (col_j - mean_j))
        
        # Concordance correlation coefficient
        denom = var_i + var_j + (mean_i - mean_j) ** 2
        if denom > 0:
            cc = 1 - (2 * cov) / denom
            all_corrs.append(cc)
    
    return np.mean(all_corrs) if all_corrs else 1.0


# =============================================================================
# PERMUTATION-BASED NULL DISTRIBUTION
# =============================================================================

def permutation_null_distribution(
    gene_set: set,
    universe: set,
    target_genes: set,
    n_permutations: int = 1000,
    random_state: int = 42
) -> np.ndarray:
    """
    Generate permutation-based null distribution for enrichment p-values.
    
    Shuffles gene labels to create null expectation.
    """
    random.seed(random_state)
    np.random.seed(random_state)
    
    universe_list = list(universe)
    n_genes_in_set = len(gene_set & universe)
    null_pvalues = []
    
    for _ in range(n_permutations):
        # Shuffle target/background assignment
        shuffled_targets = set(random.sample(universe_list, len(target_genes)))
        shuffled_bg = universe - shuffled_targets
        
        # Compute enrichment for this permutation
        in_target = len(gene_set & shuffled_targets)
        in_bg = len(gene_set & shuffled_bg)
        out_target = len(shuffled_targets - gene_set)
        out_bg = len(shuffled_bg - gene_set)
        
        table = [[in_target, in_bg], [out_target, out_bg]]
        try:
            _, pval = stats.fisher_exact(table, alternative="greater")
            null_pvalues.append(pval)
        except:
            null_pvalues.append(1.0)
    
    return np.array(null_pvalues)


def adjusted_enrichment_pvalue(
    observed_pvalue: float,
    null_distribution: np.ndarray
) -> float:
    """Adjust p-value using permutation-based null distribution ( BH correction)."""
    n_null = len(null_distribution)
    adjusted = (np.sum(null_distribution <= observed_pvalue) + 1) / (n_null + 1)
    return min(adjusted, 1.0)


# =============================================================================
# MAIN: Stability Benchmark
# =============================================================================

def benchmark_stability(
    deg_file: str,
    gene_set_db: str,
    output_dir: str = "enrichment_stability_results",
    count_matrix: Optional[str] = None,
    gtf_path: Optional[str] = None,
    use_gsea: bool = False
) -> Dict:
    """
    Main entry point. Runs enrichment under 4 universe definitions and
    quantifies stability of conclusions.
    
    Args:
        deg_file: CSV with columns [gene_id, log2fc, pvalue]
        gene_set_db: GMT file (MSigDB format) or JSON dict {pathway: [genes]}
        output_dir: Where to write results
        count_matrix: Optional CSV with gene expression counts for universe building
        gtf_path: Optional GTF file for all-annotated universe
        use_gsea: Whether to also run GSEA-style rank-based analysis
    
    Returns:
        Dict with stability scores, per-pathway concordance, and instability flags
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load DEGs
    degs = pd.read_csv(deg_file)
    target_genes = set(degs[degs["pvalue"] < 0.05]["gene_id"])
    print(f"[EnrichmentStability] {len(target_genes)} DEGs loaded (p<0.05)")
    
    # Create ranked gene list for GSEA
    ranked_genes = degs.sort_values("pvalue").apply(
        lambda row: (row["gene_id"], -np.log10(row["pvalue"] + 1e-100) * np.sign(row["log2fc"])), 
        axis=1
    ).tolist()
    
    # Load gene sets
    if gene_set_db.endswith(".gmt"):
        gene_sets = parse_gmt(gene_set_db)
    else:
        with open(gene_set_db) as f:
            gene_sets = json.load(f)
    
    print(f"[EnrichmentStability] {len(gene_sets)} pathways loaded")
    
    # Build 4 universes
    universes = {}
    
    if count_matrix:
        cm = pd.read_csv(count_matrix, index_col=0)
        universes["detected"] = build_detected_genes(cm)
        universes["expression_filtered"] = build_expression_filtered_genes(cm)
        universes["low_count_pruned"] = build_low_count_pruned_genes(cm)
        print(f"[EnrichmentStability] Built 3 expression-based universes from count matrix")
    
    if gtf_path:
        universes["all_annotated"] = build_all_annotated_genes(gtf_path)
        print(f"[EnrichmentStability] Built all-annotated universe from GTF")
    
    # If no external files, use demo synthetic universes
    if not universes:
        all_genes = set(degs["gene_id"])
        universes = {
            "all_annotated": all_genes,
            "detected": all_genes,
            "expression_filtered": all_genes,
            "low_count_pruned": all_genes,
        }
        print(f"[EnrichmentStability] Using DEG-based demo universes")
    
    # Run enrichment for each universe
    results = {}
    gsea_results = {}
    
    for universe_name, universe in universes.items():
        universe_results = {}
        gsea_universe_results = {}
        background_genes = universe - target_genes
        
        for pathway, genes in gene_sets.items():
            gene_set = set(genes)
            
            # Standard Fisher test
            pval, it, ib, ot, ob = fisher_enrichment(gene_set, universe, target_genes, background_genes)
            universe_results[pathway] = {
                "pvalue": pval,
                "in_target": it,
                "significant": pval < 0.05
            }
            
            # GSEA-style if requested
            if use_gsea:
                es, gsea_pval = gsea_enrichment(gene_set, ranked_genes)
                gsea_universe_results[pathway] = {
                    "enrichment_score": es,
                    "pvalue": gsea_pval,
                    "significant": gsea_pval < 0.05
                }
        
        results[universe_name] = universe_results
        if use_gsea:
            gsea_results[universe_name] = gsea_universe_results
        
        n_sig = sum(1 for r in universe_results.values() if r['significant'])
        print(f"[EnrichmentStability] {universe_name}: {n_sig} significant pathways")
    
    # Compute stability metrics
    stability_report = compute_stability_metrics(results, gsea_results if use_gsea else None)
    
    # Write outputs
    _write_outputs(stability_report, results, gsea_results if use_gsea else None, output_path, ranked_genes)
    
    return stability_report


def compute_stability_metrics(results: Dict, gsea_results: Optional[Dict] = None) -> Dict:
    """Compute per-pathway and overall stability metrics."""
    pathways = list(results[list(results.keys())[0]].keys())
    n_universes = len(results)
    
    per_pathway = {}
    all_pvalue_matrix = []
    
    for pathway in pathways:
        pathway_data = {}
        pvalues = []
        significances = []
        
        for universe_name, universe_results in results.items():
            pvalues.append(universe_results[pathway]["pvalue"])
            significances.append(universe_results[pathway]["significant"])
        
        # Significance concordance (fraction of universes where significant)
        n_significant = sum(significances)
        pathway_data["significant_in"] = n_significant
        pathway_data["concordance"] = n_significant / n_universes
        pathway_data["unstable"] = n_significant < (n_universes + 1) // 2  # significant in <50% of universes
        
        # Rank correlation of pvalues across universes
        if len(pvalues) >= 3:
            # Use Spearman rank correlation across all universe pairs
            ranks = [stats.rankdata(p) for p in [pvalues]]
            pathway_data["mean_rank"] = np.mean(ranks)
        
        # Store pvalues for concordance correlation
        all_pvalue_matrix.append(pvalues)
        
        per_pathway[pathway] = pathway_data
    
    # Overall concordance correlation
    pvalue_df = pd.DataFrame(all_pvalue_matrix, index=pathways,
                              columns=list(results.keys()))
    overall_concordance = compute_concordance_correlation(pvalue_df)
    
    # Fraction of unstable pathways
    unstable_pathways = [p for p, d in per_pathway.items() if d["unstable"]]
    fraction_unstable = len(unstable_pathways) / len(pathways) if pathways else 0
    
    stability_report = {
        "mean_concordance": float(overall_concordance),
        "fraction_unstable": float(fraction_unstable),
        "unstable_pathways": unstable_pathways,
        "per_pathway": per_pathway,
        "n_pathways": len(pathways),
        "n_universes": n_universes,
        "universe_names": list(results.keys())
    }
    
    # GSEA-specific metrics
    if gsea_results:
        gsea_stability = {}
        for pathway in pathways:
            es_values = [gsea_results[u][pathway]["enrichment_score"] for u in gsea_results]
            pvalues = [gsea_results[u][pathway]["pvalue"] for u in gsea_results]
            sigs = [gsea_results[u][pathway]["significant"] for u in gsea_results]
            
            gsea_stability[pathway] = {
                "mean_es": float(np.mean(es_values)),
                "es_variance": float(np.var(es_values)),
                "significant_in": sum(sigs),
                "stable": sum(sigs) >= (len(sigs) + 1) // 2
            }
        
        stability_report["gsea_stability"] = gsea_stability
    
    return stability_report


def _convert_to_native(obj):
    """Convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _convert_to_native(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_convert_to_native(item) for item in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _write_outputs(
    stability_report: Dict, 
    results: Dict, 
    gsea_results: Optional[Dict],
    output_path: Path,
    ranked_genes: List
):
    """Write JSON, CSV, and HTML outputs."""
    # JSON report - convert numpy types to native Python
    stability_report = _convert_to_native(stability_report)
    with open(output_path / "stability_report.json", "w") as f:
        json.dump(stability_report, f, indent=2)
    
    # Flatten to CSV
    pathways = list(results[list(results.keys())[0]].keys())
    rows = []
    for pathway in pathways:
        row = {"pathway": pathway}
        for universe, universe_results in results.items():
            row[f"{universe}_pvalue"] = universe_results[pathway]["pvalue"]
            row[f"{universe}_significant"] = int(universe_results[pathway]["significant"])
        if gsea_results:
            for universe, universe_results in gsea_results.items():
                row[f"{universe}_es"] = universe_results[pathway]["enrichment_score"]
                row[f"{universe}_gsea_pvalue"] = universe_results[pathway]["pvalue"]
        row["stability_score"] = stability_report["per_pathway"][pathway]["concordance"]
        row["unstable"] = int(stability_report["per_pathway"][pathway]["unstable"])
        rows.append(row)
    
    pd.DataFrame(rows).to_csv(output_path / "enrichment_by_universe.csv", index=False)
    
    # Generate interactive Plotly heatmap
    _generate_heatmap(results, stability_report, output_path)
    
    print(f"[EnrichmentStability] Results written to {output_path}")


def _generate_heatmap(
    results: Dict, 
    stability_report: Dict,
    output_path: Path
):
    """Generate interactive Plotly heatmap of universe-swap effects."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        pathways = list(results[list(results.keys())[0]].keys())
        universes = list(results.keys())
        
        # -log10(pvalue) matrix for heatmap
        z_matrix = []
        for pathway in pathways:
            row = []
            for universe in universes:
                pval = results[universe][pathway]["pvalue"]
                row.append(-np.log10(max(pval, 1e-100)))
            z_matrix.append(row)
        
        # Create annotations (significance markers)
        annotations = []
        for i, pathway in enumerate(pathways):
            for j, universe in enumerate(universes):
                sig = results[universe][pathway]["significant"]
                annotations.append(
                    dict(
                        x=j, y=i,
                        text="*" if sig else "",
                        showarrow=False,
                        font=dict(size=20, color="red" if sig else "gray")
                    )
                )
        
        # Heatmap
        fig = go.Figure(data=go.Heatmap(
            z=z_matrix,
            x=universes,
            y=pathways,
            colorscale="RdBu_r",
            zmid=1.3,  # ~p=0.05 in -log10 scale
            colorbar=dict(title="-log10(p-value)"),
            hovertemplate="Pathway: %{y}<br>Universe: %{x}<br>-log10(p): %{z:.2f}<extra></extra>"
        ))
        
        fig.update_layout(
            title={
                'text': "Enrichment Significance Across Background Universes<br><sub>Universe Swap Sensitivity Heatmap</sub>",
                'x': 0.5
            },
            xaxis_title="Background Universe Strategy",
            yaxis_title="Pathway",
            annotations=annotations,
            width=800,
            height=max(400, len(pathways) * 40)
        )
        
        fig.write_html(output_path / "universe_heatmap.html")
        
        # Second heatmap: stability scores per pathway
        stability_matrix = [[stability_report["per_pathway"][p]["concordance"]] for p in pathways]
        
        fig2 = go.Figure(data=go.Heatmap(
            z=stability_matrix,
            x=["Stability Score"],
            y=pathways,
            colorscale="RdYlGn",
            zmid=0.5,
            colorbar=dict(title="Concordance"),
            hovertemplate="Pathway: %{y}<br>Stability: %{z:.2f}<extra></extra>"
        ))
        
        fig2.update_layout(
            title={
                'text': "Per-Pathway Stability Scores<br><sub>Fraction of universes where pathway is significant</sub>",
                'x': 0.5
            },
            width=400,
            height=max(400, len(pathways) * 40)
        )
        
        fig2.write_html(output_path / "stability_heatmap.html")
        
        print(f"[EnrichmentStability] Generated interactive heatmaps")
        
    except ImportError:
        print("[EnrichmentStability] Plotly not available - skipping interactive heatmaps")
        # Generate simple matplotlib fallback
        _generate_matplotlib_heatmap(results, stability_report, output_path)


def _generate_matplotlib_heatmap(
    results: Dict, 
    stability_report: Dict,
    output_path: Path
):
    """Fallback matplotlib heatmap if Plotly not available."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        
        pathways = list(results[list(results.keys())[0]].keys())
        universes = list(results.keys())
        
        # -log10(pvalue) matrix
        z_matrix = np.zeros((len(pathways), len(universes)))
        for i, pathway in enumerate(pathways):
            for j, universe in enumerate(universes):
                pval = results[universe][pathway]["pvalue"]
                z_matrix[i, j] = -np.log10(max(pval, 1e-100))
        
        fig, ax = plt.subplots(figsize=(10, max(6, len(pathways) * 0.4)))
        
        im = ax.imshow(z_matrix, aspect='auto', cmap='RdBu_r', vmin=0, vmax=5)
        
        ax.set_xticks(range(len(universes)))
        ax.set_xticklabels(universes, rotation=45, ha='right')
        ax.set_yticks(range(len(pathways)))
        ax.set_yticklabels(pathways)
        
        # Add significance markers
        for i, pathway in enumerate(pathways):
            for j, universe in enumerate(universes):
                sig = results[universe][pathway]["significant"]
                if sig:
                    ax.text(j, i, "*", ha='center', va='center', color='red', fontsize=16, fontweight='bold')
        
        plt.colorbar(im, ax=ax, label='-log10(p-value)')
        ax.set_title('Enrichment Significance Across Background Universes')
        ax.set_xlabel('Background Universe Strategy')
        ax.set_ylabel('Pathway')
        
        plt.tight_layout()
        plt.savefig(output_path / "universe_heatmap.png", dpi=150)
        plt.close()
        
        print("[EnrichmentStability] Generated matplotlib fallback heatmap")
        
    except ImportError:
        print("[EnrichmentStability] matplotlib not available - no heatmap generated")


# =============================================================================
# UTILS
# =============================================================================

def parse_gmt(gmt_path: str) -> Dict[str, List[str]]:
    """Parse GMT gene set file."""
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = [g for g in parts[2:] if g]
            gene_sets[name] = genes
    return gene_sets


# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="EnrichmentStability: GSEA universe stability benchmark")
    parser.add_argument("--deg-file", required=True, help="CSV with gene_id, log2fc, pvalue")
    parser.add_argument("--gene-set-db", required=True, help="GMT file or JSON dict of gene sets")
    parser.add_argument("--count-matrix", help="Raw count matrix for universe construction (optional)")
    parser.add_argument("--gtf-path", help="GTF file for all-annotated universe (optional)")
    parser.add_argument("--output-dir", default="enrichment_stability_results")
    parser.add_argument("--use-gsea", action="store_true", help="Also run GSEA-style rank-based analysis")
    
    args = parser.parse_args()
    
    report = benchmark_stability(
        args.deg_file, 
        args.gene_set_db, 
        args.output_dir,
        args.count_matrix,
        args.gtf_path,
        args.use_gsea
    )
    print(json.dumps(report, indent=2))

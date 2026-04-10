#!/usr/bin/env python3
"""
Demo for EnrichmentStability.

Generates synthetic DEG data and gene set database to demonstrate the workflow.
Uses sample SARS-CoV-2 and E. coli gene annotations for biological realism.
Expected output: stability_report.json + enrichment_by_universe.csv + heatmap.html
"""

import json
import numpy as np
import pandas as pd
import sys
import os
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from EnrichmentStability import (
    benchmark_stability,
    build_detected_genes,
    build_expression_filtered_genes,
    build_low_count_pruned_genes,
    fisher_enrichment,
    parse_gmt
)


def generate_sars_cov2_demo_data():
    """
    Generate demo dataset based on SARS-CoV-2 gene annotations.
    Realistic gene sets from KEGG coronavirus pathway and host response.
    """
    np.random.seed(42)
    
    # SARS-CoV-2 and human host response genes (realistic gene IDs)
    sars_genes = [
        "ORF1a", "ORF1b", "S", "E", "M", "N",  # Viral structural
        "ORF3a", "ORF3b", "ORF3c", "ORF3d",     # ORF3
        "ORF4a", "ORF4b", "ORF5",               # ORF4/5
        "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "ORF9c", "ORF10",  # Accessory
        # Human host genes involved in viral response
        "ACE2", "TMPRSS2", "ADAM17", "BSG", "HSPA5", "XPO1",
        "DDX3X", "DDX5", "EIF3A", "EIF4A1", "EIF4G1",
        "RPS6", "RPL3", "RPL4", "RPL7", "RPL13",
        "STAT1", "STAT2", "STAT3", "JAK1", "TYK2",
        "IFNA1", "IFNB1", "IFNL1", "IFNL2",
        "IRF3", "IRF7", "IRF9", "NFKB1", "RELA",
        "IL6", "IL1B", "TNF", "CXCL8", "CCL2",
        "DDIT4", "GADD45A", "ATF4", "ATF3",
        "HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1",
        "DDX21", "DHX36", "MOV10", "LIN28A",
        "AKT1", "AKT2", "MTOR", "RPS6KB1",
        "MAPK1", "MAPK3", "MAPK14", "JUN", "FOS",
        "TP53", "MDM2", "BCL2", "BAX", "CASP3",
        "CTNNB1", "APC", "GSK3B",
        "LDHA", "LDHB", "PKM", "PGK1", "ENO1",
        "COX1", "COX2", "ATP5A1", "ATP5B",
        "SOD1", "SOD2", "CAT", "GPX1", "NFE2L2",
        "NLRP3", "AIM2", "PYCARD", "CASP1",
        "CGAS", "STING1", "TBK1", "IKBKE",
        "SUMO1", "SUMO2", "SENP1", "SENP2",
        "HDAC1", "HDAC2", "HDAC3", "SIRT1",
        "PARP1", "PARP2", "BRCA1", "BRCA2",
        "RAD51", "MRE11", "NBN", "ATM",
        "CHEK1", "CHEK2", "CDC25A", "CDC25B",
        "CDK1", "CDK2", "CDK4", "CDK6", "CCNA2", "CCNB1", "CCND1"
    ]
    
    # E. coli genes for gut microbiome demo
    ecoli_genes = [
        # Central metabolism
        "gyrA", "gyrB", "topA", "recA", "lexA", "umuC", "umuD",
        "dnaA", "dnaB", "dnaC", "dnaG", "dnaN", "dnaQ",
        "rpoA", "rpoB", "rpoC", "rpoD", "rpoH", "rpoE", "rpoS",
        "rpsA", "rpsB", "rpsC", "rpsD", "rpsE", "rpsF", "rpsG", "rpsH", "rpsI", "rpsJ", "rpsK", "rpsL", "rpsM", "rpsN", "rpsO", "rpsP", "rpsQ", "rpsR", "rpsS", "rpsT", "rpsU",
        "rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplJ", "rplK", "rplL", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU", "rplV", "rplW", "rplX", "rplY",
        "lpp", "ompA", "ompF", "ompC", "lamB", "malE", "malK",
        "galK", "galT", "galE", "galR", "galS",
        "lacI", "lacZ", "lacY", "lacA", "lacO",
        "araC", "araA", "araB", "araD", "araE", "araF",
        "trpA", "trpB", "trpC", "trpD", "trpE", "trpR",
        "hisA", "hisB", "hisC", "hisD", "hisF", "hisG", "hisH", "hisI", "hisJ", "hisK", "hisM", "hisP", "hisQ", "hisR", "hisS", "hisU", "hisW",
        "leuA", "leuB", "leuC", "leuD", "leuE", "leuL", "leuO", "leuQ", "leuZ",
        "metA", "metB", "metC", "metE", "metF", "metH", "metK", "metL", "metR",
        "argA", "argB", "argC", "argD", "argE", "argF", "argG", "argH", "argI", "argJ", "argK", "argM", "argO", "argP", "argQ", "argR", "argS", "argU", "argV", "argW", "argX", "argY", "argZ",
        "ilvA", "ilvB", "ilvC", "ilvD", "ilvE", "ilvG", "ilvH", "ilvI", "ilvK", "ilvL", "ilvM", "ilvN", "ilvO", "ilvP", "ilvQ", "ilvR", "ilvS", "ilvX", "ilvY", "ilvZ",
        "thrA", "thrB", "thrC", "thrL", "thrR",
        "adeB", "adeC", "adeK", "adeA", "adeP", "adeR",
        "gltA", "gltB", "gltD", "gltF", "gltJ", "gltK", "gltL", "gltX",
        "sucA", "sucB", "sucC", "sucD", "sucE", "sucF", "sucG", "sucH", "sucI", "sucJ", "sucK",
        "mdh", "icdA", "fumA", "fumB", "fumC", "fumD", "fumE", "fumF", "fumG", "fumH", "fumI", "fumK", "fumM",
        "eno", "pgk", "gpmA", "gpmB", "gpmI", "tpiA",
        "pykA", "pykF", "zbwA",
        "pgi", "pgl", "pgm", "galU", "galM", "galK", "galT", "galR",
        "glmS", "glmM", "glmU",
        "murA", "murB", "murC", "murD", "murE", "murF", "mraY",
        "ftsA", "ftsB", "ftsI", "ftsL", "ftsQ", "ftsW", "ftsX", "ftsY", "ftsZ",
        "lon", "clpP", "clpX", "clpA", "clpB", "hslU", "hslV",
        "groEL", "groES", "dnaK", "dnaJ", "grpE",
        "secA", "secB", "secD", "secE", "secF", "secG", "secY",
        "tatA", "tatB", "tatC", "tatD", "tatE", "tatG", "tatH",
        "lptA", "lptB", "lptC", "lptD", "lptE", "lptF", "lptG", "lptH", "lptJ",
        "degP", "degQ", "degS", "sapA", "sapB", "sapC", "sapD", "sapF",
        "surA", "pal", "Skp", "OmpLA", "OmpLB", "OmpLC",
        "msbA", "lpt", "arnA", "arnB", "arnC", "arnD", "arnE", "arnF", "arnT",
        "waaA", "waaB", "waaC", "waaD", "waaE", "waaF", "waaG", "waaH", "waaJ", "waaK", "waaL", "waaO", "waaP", "waaQ", "waaR", "waaS", "waaU", "waaV", "waaW", "waaX", "waaY", "waaZ",
        "waeA", "waeI", "waeK", "waeL", "waeP", "waeQ", "waeR", "waeW", "waeX", "waeY", "waeZ",
        "lpsA", "lpsB", "lpsC", "lpsD", "lpsE", "lpsX"
    ]
    
    # Combine all genes into universe
    all_genes = list(set(sars_genes + ecoli_genes))
    np.random.shuffle(all_genes)
    
    # Generate synthetic DEGs with realistic distribution
    n_degs = 60
    deg_genes = np.random.choice(all_genes, size=n_degs, replace=False)
    log2fc = np.concatenate([
        np.random.uniform(1, 3, size=n_degs//2),      # Upregulated
        np.random.uniform(-3, -1, size=n_degs//2)     # Downregulated
    ])
    np.random.shuffle(log2fc)
    
    # Realistic p-value distribution (many near 1, few very significant)
    pvalue = np.random.exponential(0.15, size=n_degs)
    pvalue = np.clip(pvalue, 0.0001, 1)
    
    deg_df = pd.DataFrame({
        "gene_id": deg_genes,
        "log2fc": log2fc,
        "pvalue": pvalue
    })
    deg_df.to_csv("/tmp/EnrichmentStability_demo_DEGs.csv", index=False)
    
    # Generate realistic gene sets (pathways)
    gene_sets = {
        # SARS-CoV-2 related
        "VIRAL_ENTRY": ["ACE2", "TMPRSS2", "ADAM17", "BSG", "CTNNB1", "AKT1", "MAPK1", "MAPK3"],
        "INTERFERON_RESPONSE": ["IFNA1", "IFNB1", "IFNL1", "STAT1", "STAT2", "JAK1", "TYK2", "IRF3", "IRF7", "IRF9"],
        "INFLAMMATION": ["IL6", "IL1B", "TNF", "CXCL8", "CCL2", "NFKB1", "RELA", "MAPK14", "JUN", "FOS"],
        "VIRAL_REPLICATION": ["ORF1a", "ORF1b", "S", "E", "M", "N", "ORF3a", "ORF6", "ORF7a", "ORF8", "DDX3X", "DDX5", "EIF3A"],
        "STRESS_RESPONSE": ["HSPA1A", "HSPA5", "HSP90AA1", "DDIT4", "ATF4", "ATF3", "XPO1", "SUMO1"],
        # E. coli / bacterial related
        "DNA_REPAIR": ["recA", "lexA", "umuC", "radA", "uvrA", "uvrB", "uvrC", "mutS", "mutL", "mutH", "dam", "topA"],
        "TRANSCRIPTION": ["rpoA", "rpoB", "rpoC", "rpoD", "rpoH", "rpoE", "rpoS", "nusA", "nusG"],
        "RIBOSOME_BIOGENESIS": ["rpsA", "rpsB", "rpsC", "rpsD", "rpsE", "rplA", "rplB", "rplC", "rplD", "rplE", "rplJ", "rplK", "rplL"],
        "CENTRAL_METABOLISM": ["gltA", "sucA", "sucB", "mdh", "icdA", "fumA", "eno", "pgk", "pykA", "pgi"],
        "OUTER_MEMBRANE_BIOSYNTHESIS": ["lptA", "lptB", "lptC", "lptD", "lptE", "lptF", "lptG", "lpsA", "lpsB", "murA", "murB", "murC", "murD"],
    }
    
    with open("/tmp/EnrichmentStability_demo_genesets.json", "w") as f:
        json.dump(gene_sets, f)
    
    # Generate synthetic count matrix for universe building demos
    n_samples = 6
    n_genes = len(all_genes)
    count_data = np.random.negative_binomial(5, 0.5, size=(n_genes, n_samples))
    # Make some genes very low expressed
    low_expr_idx = np.random.choice(n_genes, size=n_genes//4, replace=False)
    for idx in low_expr_idx:
        count_data[idx, :] = np.random.poisson(0.5, size=n_samples)
    
    count_df = pd.DataFrame(count_data, index=all_genes, columns=[f"Sample_{i+1}" for i in range(n_samples)])
    count_df.to_csv("/tmp/EnrichmentStability_demo_counts.csv")
    
    print("[Demo] Generated realistic SARS-CoV-2 + E. coli dataset:")
    print(f"  - /tmp/EnrichmentStability_demo_DEGs.csv ({len(deg_df)} genes, p<0.05)")
    print(f"  - /tmp/EnrichmentStability_demo_genesets.json ({len(gene_sets)} pathways)")
    print(f"  - /tmp/EnrichmentStability_demo_counts.csv ({n_genes} genes x {n_samples} samples)")
    
    return "/tmp/EnrichmentStability_demo_DEGs.csv", "/tmp/EnrichmentStability_demo_genesets.json", "/tmp/EnrichmentStability_demo_counts.csv"


def test_universe_builders(count_matrix_path):
    """Test the 4 universe construction strategies."""
    print("\n[Demo] Testing universe construction strategies:")
    cm = pd.read_csv(count_matrix_path, index_col=0)
    
    detected = build_detected_genes(cm, detection_threshold=1.0)
    print(f"  - Detected genes (threshold=1): {len(detected)}")
    
    expr_filtered = build_expression_filtered_genes(cm, min_samples=3, min_counts=10)
    print(f"  - Expression-filtered (3 samples, 10 counts): {len(expr_filtered)}")
    
    low_pruned = build_low_count_pruned_genes(cm, percentile_cutoff=25)
    print(f"  - Low-count pruned (bottom 25%): {len(low_pruned)}")
    
    return {
        "detected": detected,
        "expression_filtered": expr_filtered,
        "low_count_pruned": low_pruned,
        "all_annotated": set(cm.index)
    }


def test_enrichment_analysis(deg_file, gene_set_file, count_matrix_path):
    """Run full enrichment stability benchmark."""
    print("\n[Demo] Running enrichment stability benchmark...")
    
    # Run benchmark
    report = benchmark_stability(
        deg_file=deg_file,
        gene_set_db=gene_set_file,
        output_dir="/tmp/EnrichmentStability_results",
        count_matrix=count_matrix_path,
        use_gsea=True
    )
    
    return report


def main():
    print("=" * 70)
    print("EnrichmentStability Demo - SARS-CoV-2 / E. coli Analysis")
    print("=" * 70)
    
    # Generate demo data
    deg_file, gs_file, count_file = generate_sars_cov2_demo_data()
    
    # Test universe builders
    universes = test_universe_builders(count_file)
    
    # Run analysis
    report = test_enrichment_analysis(deg_file, gs_file, count_file)
    
    print("\n" + "=" * 70)
    print("STABILITY REPORT SUMMARY")
    print("=" * 70)
    print(f"\nMean Concordance: {report['mean_concordance']:.3f}")
    print(f"Unstable Pathways: {len(report['unstable_pathways'])} / {report['n_pathways']}")
    print(f"Fraction Unstable: {report['fraction_unstable']:.1%}")
    
    if report['unstable_pathways']:
        print(f"\nUnstable pathways (significant in <50% of universes):")
        for p in report['unstable_pathways']:
            print(f"  - {p}")
    
    print("\nPer-Pathway Stability:")
    for pathway, data in report['per_pathway'].items():
        status = "UNSTABLE" if data['unstable'] else "stable"
        print(f"  {pathway}: {data['significant_in']}/{data['concordance']:.0%} ({status})")
    
    print("\n" + "=" * 70)
    print("EXPECTED OUTPUT FILES:")
    print("=" * 70)
    print("  /tmp/EnrichmentStability_results/")
    print("  ├── stability_report.json      # Full stability metrics")
    print("  ├── enrichment_by_universe.csv # Per-universe p-values")
    print("  ├── universe_heatmap.html       # Interactive Plotly heatmap")
    print("  └── stability_heatmap.html     # Stability score heatmap")
    
    print("\n[Demo] Complete!")
    return report


if __name__ == "__main__":
    main()

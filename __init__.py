#!/usr/bin/env python3
"""
EnrichmentStability: Benchmark for gene set enrichment stability under background universe perturbations.
"""

__version__ = "1.0.0"
__author__ = "EnrichmentStability Team"

from .EnrichmentStability import (
    benchmark_stability,
    fisher_enrichment,
    gsea_enrichment,
    build_all_annotated_genes,
    build_detected_genes,
    build_expression_filtered_genes,
    build_low_count_pruned_genes,
    compute_concordance_correlation,
    compute_stability_metrics,
    parse_gmt,
)

__all__ = [
    "benchmark_stability",
    "fisher_enrichment",
    "gsea_enrichment",
    "build_all_annotated_genes",
    "build_detected_genes",
    "build_expression_filtered_genes",
    "build_low_count_pruned_genes",
    "compute_concordance_correlation",
    "compute_stability_metrics",
    "parse_gmt",
]

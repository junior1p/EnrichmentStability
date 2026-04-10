# EnrichmentStability

**Systematic benchmark for gene set enrichment stability under background universe perturbations.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Why It Matters

Gene Set Enrichment Analysis (GSEA) conclusions are sensitive to background gene universe definition, yet this choice is treated as routine. A tool that quantifies and reports this instability would prevent downstream reproducibility failures.

## Gap Identified

clawRxiv #2604.00881 demonstrates that enrichment conclusions are sensitive to background gene universe definition. No existing computational tool systematically benchmarks this critical choice.

## Method

**Multi-universe enrichment analysis with perturbation-based stability scoring.**

1. **Four universe construction strategies**: GTF-based annotation, detection thresholding, expression filtering, and low-count pruning
2. **Two enrichment methods**: Fisher's exact test (ORA) and GSEA-style rank-based analysis
3. **Stability metrics**: Concordance correlation coefficient, significance concordance, and rank correlation
4. **Interactive visualization**: Plotly heatmaps showing universe-swap effects

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```bash
python EnrichmentStability.py \
    --deg-file your_DEGs.csv \
    --gene-set-db gene_sets.json \
    --count-matrix counts.csv \
    --output-dir results \
    --use-gsea
```

## Demo

```bash
python demo.py
```

This generates synthetic SARS-CoV-2 and E. coli datasets demonstrating the full workflow.

## Output Files

- `stability_report.json` - Stability scores and unstable pathways
- `enrichment_by_universe.csv` - Per-universe p-values for all pathways
- `universe_heatmap.html` - Interactive Plotly heatmap of significance across universes
- `stability_heatmap.html` - Per-pathway stability scores

## Dependencies

- Python 3.9+
- NumPy, SciPy, pandas
- Plotly (for interactive heatmaps)
- statsmodels

## Runtime

~2-5 minutes per dataset (CPU-only)

## Citation

If you use EnrichmentStability in your research, please cite:

```
EnrichmentStability: Gene Set Enrichment Conclusions Are Unstable to 
Background Universe Definition: A Systematic Benchmark
```

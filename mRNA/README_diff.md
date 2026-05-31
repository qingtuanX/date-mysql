# mRNA Differential Expression Analysis

## Overview

Differential expression analysis for three GEO datasets related to tuberculosis/lung disease using **limma-voom** pipeline.

## Datasets

| Dataset | Samples | Groups | Method | Input |
|---------|---------|--------|--------|-------|
| **GSE148036** | 20 | Normal (5), Tuberculosis (5), Adenocarcinoma (5), Sacrodosis (5) | limma-voom | Raw counts |
| **GSE158767** | 10 | H (5), L (5) | limma-trend | FPKM |
| **GSE99374** | 40 | Healthy (21), Latent TB infection (19) | limma-voom | Raw counts |

## Pipeline

```
Raw counts → CPM filter → TMM normalization → voom transformation
         ↓
    limma lmFit → contrasts.fit (Target vs Others) → eBayes
         ↓
  topTable → log2FC, PValue, FDR sorted by B-statistic
```

### Key Features

- **GSE148036 / GSE99374**: `edgeR::DGEList` → `calcNormFactors` → `limma::voom()` → `lmFit` → `eBayes`
- **GSE158767**: FPKM → log2(FPKM + 0.5) → `lmFit` → `eBayes(trend = TRUE)` (limma-trend for normalized data)
- Low-expression filtering: CPM > 1 in ≥ 2 samples (GSE148036/99374), mean FPKM ≥ 1 (GSE158767)
- Multiple testing correction: Benjamini-Hochberg FDR

## Output Files

```
差异化分析/
├── GSE148036/
│   ├── GSE148036_limma_test_result.*_vs_Others.txt   # Full results (tab-separated)
│   ├── plots/                                         # Volcano + Heatmap PNG
│   └── run_limma_gse148036.R                          # Analysis script
├── GSE158767/
│   ├── GSE158767_limma_test_result.*_vs_Others.txt
│   ├── plots/
│   └── run_limma_gse158767.R
└── GSE99374/
    ├── GSE99374_limma_test_result.*_vs_Others.txt
    ├── plots/
    └── run_limma_gse99374.R
```

### Result File Columns

| Column | Description |
|--------|-------------|
| `log2FC` | Log2 fold change (Target vs Others) |
| `t` | Moderated t-statistic |
| `B` | B-statistic (log-odds of differential expression) |
| `PValue` | Raw P-value |
| `FDR` | Benjamini-Hochberg adjusted P-value |

## Visualizations

- **Volcano plot**: log2FC vs -log10(PValue), with |log2FC| ≥ 1 & P < 0.05 highlighted
- **Heatmap**: Top 50 DEGs by PValue, z-score scaled, Ward.D2 clustering

## Notes

- Contrast design is **"one group vs all others"** (not pairwise). This provides a broader baseline comparison.
- For GSE148036, pairwise comparisons (e.g., AD vs Normal only) are available in `../result/result/gse148036_results/`.
- CPM filter removes genes with consistently very low expression to improve statistical power.

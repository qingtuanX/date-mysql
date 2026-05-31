# TME Immune Infiltration Analysis

## Overview

Tumor microenvironment (TME) immune cell deconvolution for three GEO datasets using **IOBR** package with multiple algorithms.

## Pipeline

```
Raw counts/FPKM → TPM normalization → Gene symbol aggregation
         ↓
   IOBR::deconvo_tme (multiple methods) → Cell fraction matrices
         ↓
   Merge all methods → tme_combine.csv + heatmap
```

### TPM Conversion

- **GSE148036**: ENSG counts → get gene lengths (org.Hs.eg.db) → TPM = RPK / sum(RPK) × 10⁶
- **GSE158767**: FPKM → TPM ≈ FPKM / sum(FPKM) × 10⁶ (same species)
- **GSE99374**: Gene symbols → ENSEMBL → gene lengths → TPM

TPM is required by deconvolution methods (CIBERSORT, EPIC, quanTIseq, TIMER) for accurate cell proportion estimation.

## Methods

| Method | Description | GSE148036 | GSE158767 | GSE99374 |
|--------|-------------|:---:|:---:|:---:|
| **MCPcounter** | Marker gene-based absolute scores | ✓ | ✓ | ✓ |
| **ESTIMATE** | Stromal/Immune/ESTIMATE composite scores | ✓ | ✓ | ✓ |
| **EPIC** | Immune + stromal cell fractions (ref: RNA-seq) | ✓ | ✓ | ✓ |
| **IPS** | Immunophenoscore (MHC, EC, SC, CP) | ✓ | ✓ | ✓ |
| **CIBERSORT** | 22 immune cell types (ref: LM22 signature) | ✓ | ✓ | ✓ |
| **TIMER** | 6 immune cell types (cancer-type specific) | ✓ | ✓ | — |
| **quanTIseq** | 10 cell types, absolute fractions (lsei) | ✓ | — | — |

- **TIMER** uses `luad` (lung adenocarcinoma) reference for lung datasets; not applicable for PBMC.
- **quanTIseq** removed for GSE158767 (57% signature gene coverage), not applicable for PBMC.
- **xCell** removed from all (64 cell types, many tissue-irrelevant for lung/PBMC).

## Output Files

```
免疫浸润/
├── GSE148036/  (7 methods)
│   ├── *.csv                    # Individual method results
│   ├── tme_combine.csv          # Merged all methods
│   ├── TPM_matrix.csv           # TPM expression matrix used
│   └── plots/immune_heatmap.png # Hierarchical clustering heatmap
├── GSE158767/  (6 methods)
└── GSE99374/   (5 methods)
```

### Result File Columns (varies by method)

| Method | Key Columns |
|--------|-------------|
| CIBERSORT | B_cells_naive/memory, CD4/CD8_T_cells, NK, Macrophage M0/M1/M2, Neutrophils, etc. (22 types) |
| EPIC | Bcells, CAFs, CD4/CD8_Tcells, Endothelial, Macrophages, NKcells |
| MCPcounter | T_cells, CD8_T_cells, B_lineage, NK_cells, Monocytic_lineage, Neutrophils, Fibroblasts |
| ESTIMATE | StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity |
| TIMER | B_cell, T_cell_CD4/CD8, Neutrophil, Macrophage, DC |
| quanTIseq | B_cells, Macrophages_M1/M2, NK_cells, T_cells_CD4/CD8, Tregs, Dendritic_cells |
| IPS | MHC, EC, SC, CP, AZ, IPS |

## Interpretation Notes

- **CIBERSORT/EPIC/quanTIseq**: Cell fractions sum to ~1 (relative proportions)
- **MCPcounter**: Arbitrary units (relative abundance scores, not proportions)
- **ESTIMATE**: Composite scores (higher = more stromal/immune infiltration)
- **IPS**: Z-scores, higher = more immunogenic
- **TIMER**: Relative fractions (not constrained to sum=1)

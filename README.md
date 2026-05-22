# date-mysql: mRNA Transcriptomics Analysis Project

This project focuses on the differential expression analysis (DEA) of mRNA datasets from the Gene Expression Omnibus (GEO). It utilizes the `limma` package in R to process raw expression data and generates scientific visualizations such as Volcano Plots.

## 📁 Project Structure

- **`/mRNA`**: Core data and analysis directory.
  - **`GSE148036`, `GSE158767`, `GSE99374`**: Datasets retrieved from GEO, including metadata and expression matrices.
  - **`【378】limma进行多组差异分析`**: Contains R Markdown scripts (`Multiclasslimma.Rmd`) for performing differential analysis between multiple groups (e.g., MITF-low vs Others).
  - **`【199】多组火山图` / `【206】火山图`**: Scripts and templates for generating SCI-standard Volcano plots.
  - **`date-mysql.code-workspace`**: VS Code configuration for the project.

## 🚀 Analysis Workflow

1. **Data Acquisition**: Download GSE series from NCBI GEO (Metadata + FPKM/TPM matrices).
2. **Preprocessing**: Normalization and log-transformation of expression values.
3. **Differential Analysis**: 
   - Uses `limma` to calculate Fold Change and p-values.
   - Comparison groups identified: `immune_vs_Others`, `keratin_vs_Others`, `MITF-low_vs_Others`.
4. **Visualization**: Mapping results to Volcano plots where:
   - **X-axis**: Log2 Fold Change (Effect size).
   - **Y-axis**: -Log10 adj. P-value (Significance).

## 🛠 Prerequisites

- **R Environment**: R >= 4.0
- **Packages**: `limma`, `ggplot2`, `dplyr`, `Biobase`.
- **Database**: MySQL (intended for storing dynamic mapping coordinates).

---
*Note: This project is part of an SCI research plotting pipeline.*

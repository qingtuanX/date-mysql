# ================================================================
# GSE148036 免疫浸润热力图
# ================================================================

rm(list = ls())
setwd('d:/CODE/date-mysql/mRNA/GSE148036/')
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")

if (!file.exists("tme_combine.csv")) {
  stop("未找到 tme_combine.csv，请先运行 run_immune_infiltration_gse148036.R")
}

dat <- read.csv("tme_combine.csv", check.names = FALSE)
rownames(dat) <- dat$ID
mat <- as.matrix(dat[, -1, drop = FALSE])
mode(mat) <- "numeric"
mat <- t(mat)

pdf("免疫浸润热力图.pdf", width = 12, height = 10)
pheatmap(
  mat,
  main = "GSE148036 Immune Infiltration Heatmap",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  fontsize = 9,
  fontsize_row = 7,
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE
)
dev.off()

cat("已生成: 免疫浸润热力图.pdf\n")

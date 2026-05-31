# ================================================================
# GSE158767 免疫浸润热力图
# ================================================================

rm(list = ls())
setwd('d:/CODE/date-mysql/mRNA/GSE158767/')
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")

if (!file.exists("tme_combine.csv")) {
  stop("未找到 tme_combine.csv，请先运行 run_immune_infiltration_gse158767.R")
}

dat <- read.csv("tme_combine.csv", check.names = FALSE)
rownames(dat) <- dat$ID
mat <- as.matrix(dat[, -1, drop = FALSE])
mode(mat) <- "numeric"
mat <- t(mat)

# 清理异常值，避免hclust报NA/NaN/Inf
mat[!is.finite(mat)] <- NA
mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
if (anyNA(mat)) {
  for (i in seq_len(nrow(mat))) {
    r <- mat[i, ]
    if (anyNA(r)) {
      fill <- median(r, na.rm = TRUE)
      if (!is.finite(fill)) fill <- 0
      r[is.na(r)] <- fill
      mat[i, ] <- r
    }
  }
}

# 去掉标准化后会产生NaN的常量行
mat <- mat[apply(mat, 1, function(x) sd(x) > 0), , drop = FALSE]
if (nrow(mat) < 2 || ncol(mat) < 2) {
  stop("热力图数据不足（清理后行或列小于2）")
}

pdf("免疫浸润热力图.pdf", width = 12, height = 10)
pheatmap(
  mat,
  main = "GSE158767 Immune Infiltration Heatmap",
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

suppressPackageStartupMessages({
  library(limma)
})

options(stringsAsFactors = FALSE)

expr_file <- "GSE148036_exp.csv"
cli_file <- "GSE148036_cli.csv"
out_dir <- "."
prefix <- "GSE148036"
overwrite <- TRUE

expr <- read.csv(
  expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

cli <- read.csv(
  cli_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

sample_col_candidates <- c("source_name_ch1", "title", "geo_accession")
group_col_candidates <- c("disease:ch1", "characteristics_ch1.3")

sample_col <- sample_col_candidates[sample_col_candidates %in% colnames(cli)][1]
group_col <- group_col_candidates[group_col_candidates %in% colnames(cli)][1]

if (is.na(sample_col) || is.na(group_col)) {
  stop("临床文件缺少样本列(source_name_ch1/title/geo_accession)或分组列(disease:ch1)。")
}

subt <- data.frame(
  Sample = cli[[sample_col]],
  Group = cli[[group_col]],
  stringsAsFactors = FALSE
)

subt$Group <- sub("^disease:\\s*", "", subt$Group)
subt <- subt[!is.na(subt$Sample) & !is.na(subt$Group), , drop = FALSE]
subt <- subt[subt$Sample %in% colnames(expr), , drop = FALSE]

if (nrow(subt) == 0) {
  stop("没有可用于分析的样本匹配到表达矩阵列名。")
}

expr <- expr[, subt$Sample, drop = FALSE]
expr <- as.matrix(expr)
mode(expr) <- "numeric"

var_vec <- apply(expr, 1, sd)
threshold <- as.numeric(quantile(var_vec, probs = 0.75, na.rm = TRUE))
expr <- expr[var_vec > threshold, , drop = FALSE]

if (nrow(expr) == 0) {
  stop("方差筛选后无基因剩余，请检查输入数据。")
}

run_one_vs_others <- function(expr_mat, groups, target_group, output_file) {
  grp <- ifelse(groups == target_group, "treatment", "control")
  design <- model.matrix(~ 0 + factor(grp, levels = c("treatment", "control")))
  colnames(design) <- c("treatment", "control")

  gset <- log2(expr_mat + 1)

  fit <- lmFit(gset, design = design)
  contrasts_matrix <- makeContrasts(treatment - control, levels = design)
  fit2 <- contrasts.fit(fit, contrasts = contrasts_matrix)
  fit2 <- eBayes(fit2, 0.01)

  tt <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
  tt <- tt[, c("logFC", "t", "B", "P.Value", "adj.P.Val"), drop = FALSE]
  colnames(tt) <- c("log2FC", "t", "B", "PValue", "FDR")
  tt <- tt[order(tt$FDR), , drop = FALSE]

  write.table(tt, file = output_file, row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

groups <- subt$Group
names(groups) <- subt$Sample
all_groups <- sort(unique(groups))

cat("样本数:", ncol(expr), "\n")
cat("分组:", paste(all_groups, collapse = " | "), "\n")

for (g in all_groups) {
  out_file <- file.path(out_dir, paste0(prefix, "_limma_test_result.", safe_name(g), "_vs_Others.txt"))
  if (file.exists(out_file) && !overwrite) {
    cat("跳过已存在文件:", out_file, "\n")
    next
  }

  run_one_vs_others(expr_mat = expr, groups = groups[colnames(expr)], target_group = g, output_file = out_file)
  cat("已输出:", out_file, "\n")
}

cat("差异分析完成。\n")

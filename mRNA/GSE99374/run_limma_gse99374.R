suppressPackageStartupMessages({
  library(limma)
})

options(stringsAsFactors = FALSE)

# 输入文件
expr_file <- "G.txt"
cli_file <- "GSE99374_cli.csv"
out_dir <- "."
prefix <- "GSE99374"
overwrite <- TRUE

# 读取表达矩阵（行为基因，列为样本）
expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  row.names = 1
)

# 读取分组信息（title 对应表达矩阵样本名）
cli <- read.csv(
  cli_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

if (!("title" %in% colnames(cli)) || !("disease group:ch1" %in% colnames(cli))) {
  stop("临床文件缺少 title 或 disease group:ch1 列。")
}

subt <- data.frame(
  Sample = cli[["title"]],
  Group = cli[["disease group:ch1"]],
  stringsAsFactors = FALSE
)

subt <- subt[!is.na(subt$Sample) & !is.na(subt$Group), , drop = FALSE]
subt <- subt[subt$Sample %in% colnames(expr), , drop = FALSE]

if (nrow(subt) == 0) {
  stop("没有可用于分析的样本匹配到表达矩阵列名。")
}

# 仅保留匹配到分组信息的样本
expr <- expr[, subt$Sample, drop = FALSE]

# 为了加速，保留方差前 25% 的基因（与模板一致）
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

  # 模板使用 log2(count + 1)
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

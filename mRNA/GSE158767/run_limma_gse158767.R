suppressPackageStartupMessages({
  library(limma)
})

options(stringsAsFactors = FALSE)

expr_file <- "merged_data_dplyr.csv"
out_dir <- "."
prefix <- "GSE158767"
overwrite <- TRUE

raw_df <- read.csv(
  expr_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_cols <- grep("^[HL][0-9]+$", colnames(raw_df), value = TRUE)
if (length(sample_cols) < 2) {
  stop("未识别到样本列（如 H1/L1），请检查输入文件。")
}

if ("hgnc_symbol" %in% colnames(raw_df)) {
  gene_id <- raw_df$hgnc_symbol
  gene_id[is.na(gene_id) | gene_id == ""] <- raw_df$Transcript_id[is.na(gene_id) | gene_id == ""]
} else {
  gene_id <- raw_df$Transcript_id
}

expr_mat <- as.matrix(raw_df[, sample_cols, drop = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_id

expr <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)

subt <- data.frame(
  Sample = sample_cols,
  Group = ifelse(grepl("^H", sample_cols), "H", "L"),
  stringsAsFactors = FALSE
)

expr <- expr[, subt$Sample, drop = FALSE]

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

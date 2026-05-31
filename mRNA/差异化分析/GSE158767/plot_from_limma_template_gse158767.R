suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)

expr_file <- "merged_data_dplyr.csv"
res_pattern <- "^GSE158767_limma_test_result\\..+_vs_Others\\.txt$"
out_dir <- "plots"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
expr_log2 <- log2(expr + 1)

res_files <- list.files(path = ".", pattern = res_pattern, full.names = TRUE)
if (length(res_files) == 0) {
  stop("未找到差异结果文件，请先运行 limma 分析。")
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

extract_comp <- function(path) {
  base <- basename(path)
  sub("^GSE158767_limma_test_result\\.(.+)\\.txt$", "\\1", base)
}

draw_volcano <- function(deg, comp_name, out_png) {
  deg$gene <- rownames(deg)
  deg$neglog10p <- -log10(pmax(deg$PValue, .Machine$double.xmin))

  deg$group <- ifelse(
    deg$PValue < 0.05 & abs(deg$log2FC) >= 1,
    ifelse(deg$log2FC >= 1, "up", "down"),
    "NS"
  )
  deg$group <- factor(deg$group, levels = c("down", "NS", "up"))

  top_idx <- order(deg$PValue, decreasing = FALSE)
  top_n <- min(12, nrow(deg))
  label_genes <- deg$gene[top_idx[seq_len(top_n)]]
  deg$label <- ifelse(deg$gene %in% label_genes, deg$gene, "")

  p <- ggplot(deg, aes(x = log2FC, y = neglog10p, fill = group)) +
    geom_point(color = "black", alpha = 0.65, size = 2.5, shape = 21) +
    geom_vline(xintercept = c(-1, 1), linetype = 3, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 3, color = "black", linewidth = 0.6) +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 1000,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "black",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c(down = "#2b6cb0", NS = "#9aa0a6", up = "#d94841")) +
    labs(
      title = paste0("Volcano: ", comp_name),
      x = "log2 Fold Change",
      y = "-log10(PValue)",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "#333c41", size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "top"
    )

  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
}

draw_heatmap <- function(deg, comp_name, out_png) {
  target_group <- sub("_vs_Others$", "", comp_name)

  sample_group <- ifelse(subt$Group == target_group, target_group, "Others")
  names(sample_group) <- subt$Sample

  deg <- deg[order(deg$PValue, decreasing = FALSE), , drop = FALSE]
  genes <- intersect(rownames(deg), rownames(expr_log2))
  top_genes <- head(genes, 50)

  if (length(top_genes) < 2) {
    warning(paste0("比较 ", comp_name, " 可用于热图的基因不足，跳过热图。"))
    return(invisible(NULL))
  }

  mat <- expr_log2[top_genes, names(sample_group), drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  ann_col <- data.frame(Group = sample_group[colnames(mat)], row.names = colnames(mat))

  pheatmap(
    mat,
    filename = out_png,
    width = 9,
    height = 10,
    color = colorRampPalette(c("#2b6cb0", "white", "#d94841"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    annotation_col = ann_col,
    main = paste0("Heatmap: ", comp_name, " (Top 50 by PValue)")
  )
}

for (f in res_files) {
  comp <- extract_comp(f)
  deg <- read.table(f, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

  need_cols <- c("log2FC", "PValue", "FDR")
  if (!all(need_cols %in% colnames(deg))) {
    warning(paste0("文件缺少必要列，跳过: ", basename(f)))
    next
  }

  tag <- safe_name(comp)
  volcano_png <- file.path(out_dir, paste0("GSE158767_volcano_", tag, ".png"))
  heatmap_png <- file.path(out_dir, paste0("GSE158767_heatmap_", tag, ".png"))

  draw_volcano(deg, comp, volcano_png)
  draw_heatmap(deg, comp, heatmap_png)
  cat("已输出: ", volcano_png, "\n", sep = "")
  cat("已输出: ", heatmap_png, "\n", sep = "")
}

cat("全部作图完成。\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)

expr_file <- "G.txt"
cli_file <- "GSE99374_cli.csv"
res_pattern <- "^GSE99374_limma_test_result\\..+_vs_Others\\.txt$"
out_dir <- "plots"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read expression matrix and keep sample columns consistent with clinical table.
expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  row.names = 1
)

cli <- read.csv(
  cli_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

if (!all(c("title", "disease group:ch1") %in% colnames(cli))) {
  stop("临床文件缺少 title 或 disease group:ch1 列。")
}

subt <- data.frame(
  Sample = cli[["title"]],
  Group = cli[["disease group:ch1"]],
  stringsAsFactors = FALSE
)
subt <- subt[!is.na(subt$Sample) & !is.na(subt$Group), , drop = FALSE]
subt <- subt[subt$Sample %in% colnames(expr), , drop = FALSE]

if (nrow(subt) < 2) {
  stop("可匹配样本过少，无法作图。")
}

expr <- expr[, subt$Sample, drop = FALSE]

# Keep analysis-space genes consistent with limma script (top variance 25%).
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
  sub("^GSE99374_limma_test_result\\.(.+)\\.txt$", "\\1", base)
}

draw_volcano <- function(deg, comp_name, out_png) {
  deg$gene <- rownames(deg)
  deg$neglog10p <- -log10(pmax(deg$PValue, .Machine$double.xmin))

  # Template-style thresholds.
  deg$group <- ifelse(
    deg$PValue < 0.05 & abs(deg$log2FC) >= 1,
    ifelse(deg$log2FC >= 1, "up", "down"),
    "NS"
  )
  deg$group <- factor(deg$group, levels = c("down", "NS", "up"))

  # Always label top hits to avoid empty labels when FDR is conservative.
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

  # Recover original group label from sanitized file names.
  uniq_groups <- unique(subt$Group)
  hit <- uniq_groups[safe_name(uniq_groups) == target_group]
  if (length(hit) == 0) {
    hit <- uniq_groups[uniq_groups == gsub("_", " ", target_group)]
  }
  if (length(hit) == 0) {
    hit <- target_group
  } else {
    hit <- hit[1]
  }

  sample_group <- ifelse(subt$Group == hit, hit, "Others")
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
  volcano_png <- file.path(out_dir, paste0("GSE99374_volcano_", tag, ".png"))
  heatmap_png <- file.path(out_dir, paste0("GSE99374_heatmap_", tag, ".png"))

  draw_volcano(deg, comp, volcano_png)
  draw_heatmap(deg, comp, heatmap_png)
  cat("已输出: ", volcano_png, "\n", sep = "")
  cat("已输出: ", heatmap_png, "\n", sep = "")
}

cat("全部作图完成。\n")
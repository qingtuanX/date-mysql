suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)

# ===== Helper functions =====
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

draw_volcano <- function(deg_file, comp_name, out_png) {
  deg <- read.table(deg_file, sep = "\t", header = TRUE, row.names = 1)
  deg$gene <- rownames(deg)
  deg$neglog10p <- -log10(pmax(deg$PValue, .Machine$double.xmin))

  deg$group <- ifelse(
    deg$PValue < 0.05 & abs(deg$log2FC) >= 1,
    ifelse(deg$log2FC >= 1, "Up", "Down"), "NS"
  )
  deg$group <- factor(deg$group, levels = c("Down", "NS", "Up"))

  top_idx <- order(deg$PValue, decreasing = FALSE)
  top_n <- min(12, nrow(deg))
  deg$label <- ifelse(deg$gene %in% deg$gene[top_idx[seq_len(top_n)]], deg$gene, "")

  up_n   <- sum(deg$group == "Up")
  down_n <- sum(deg$group == "Down")

  p <- ggplot(deg, aes(x = log2FC, y = neglog10p, fill = group)) +
    geom_point(color = "black", alpha = 0.65, size = 2.5, shape = 21) +
    geom_vline(xintercept = c(-1, 1), linetype = 3, linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 3, linewidth = 0.6) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 1000,
                    box.padding = 0.5, point.padding = 0.3,
                    segment.color = "black", show.legend = FALSE) +
    scale_fill_manual(values = c(Down = "#2b6cb0", NS = "#9aa0a6", Up = "#d94841")) +
    labs(title = comp_name,
         subtitle = sprintf("Up: %d  |  Down: %d  |  NS: %d", up_n, down_n, nrow(deg)-up_n-down_n),
         x = "log2 Fold Change", y = "-log10(PValue)", fill = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "#666"),
          legend.position = "top")
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  cat(sprintf("  Volcano: %s\n", out_png))
}

draw_heatmap <- function(expr_mat, sample_group, deg_file, title, out_png, n_genes = 50) {
  deg <- read.table(deg_file, sep = "\t", header = TRUE, row.names = 1)
  deg <- deg[order(deg$PValue, decreasing = FALSE), , drop = FALSE]
  genes <- intersect(rownames(deg), rownames(expr_mat))
  top_genes <- head(genes, n_genes)
  if (length(top_genes) < 3) {
    warning(sprintf("Not enough genes for heatmap: %s", title))
    return(invisible(NULL))
  }

  mat <- expr_mat[top_genes, names(sample_group), drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  ann_col <- data.frame(Group = sample_group[colnames(mat)], row.names = colnames(mat))

  pheatmap(mat, filename = out_png, width = max(9, ncol(mat)*0.4), height = 10,
           color = colorRampPalette(c("#2b6cb0", "white", "#d94841"))(100),
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_colnames = FALSE, show_rownames = TRUE, fontsize_row = 7,
           annotation_col = ann_col, main = title)
  cat(sprintf("  Heatmap: %s\n", out_png))
}

# ========================================================
# GSE148036 plots
# ========================================================
cat("=== GSE148036 ===\n")
gse148036_dir  <- "E:/code-program/date-mysql/mRNA"
diff148036_dir <- file.path(gse148036_dir, "差异化分析/GSE148036")
data148036_dir <- file.path(gse148036_dir, "GSE148036")

# Read clinical + expression for heatmap
cli <- read.csv(file.path(data148036_dir, "GSE148036_cli.csv"), check.names = FALSE, row.names = 1)
expr <- read.csv(file.path(data148036_dir, "GSE148036_exp.csv"), check.names = FALSE, row.names = 1)
subt <- data.frame(Sample = cli[["source_name_ch1"]],
                   Group = gsub("^disease:\\s*", "", cli[["disease:ch1"]]),
                   stringsAsFactors = FALSE)
subt <- subt[subt$Sample %in% colnames(expr), ]
expr <- as.matrix(expr[, subt$Sample, drop = FALSE]); mode(expr) <- "numeric"
expr_log2 <- log2(expr + 1)

grp_map <- c("Adenocarcinoma", "Normal", "Sacrodosis", "Tuberculosis")
for (g in grp_map) {
  sname <- safe_name(g)
  f <- file.path(diff148036_dir, paste0("GSE148036_limma_test_result.", sname, "_vs_Others.txt"))
  vol_out <- file.path(diff148036_dir, "plots", paste0("GSE148036_volcano_", sname, ".png"))
  hm_out  <- file.path(diff148036_dir, "plots", paste0("GSE148036_heatmap_", sname, ".png"))

  if (file.exists(f)) {
    draw_volcano(f, paste0(g, " vs Others"), vol_out)
    sg <- ifelse(subt$Group == g, g, "Others"); names(sg) <- subt$Sample
    draw_heatmap(expr_log2, sg, f, paste0("Heatmap: ", g, " vs Others (Top 50)"), hm_out)
  }
}

# ========================================================
# GSE158767 plots
# ========================================================
cat("\n=== GSE158767 ===\n")
diff158767_dir <- file.path(gse148036_dir, "差异化分析/GSE158767")
data158767_dir <- file.path(gse148036_dir, "GSE158767")

raw <- read.csv(file.path(data158767_dir, "merged_data_dplyr.csv"), check.names = FALSE)
samp_cols <- grep("^[HL][0-9]+$", colnames(raw), value = TRUE)
raw <- raw[!is.na(raw$hgnc_symbol) & raw$hgnc_symbol != "", ]
raw$mean_val <- rowMeans(raw[, samp_cols, drop = FALSE])
best <- do.call(rbind, by(raw, raw$hgnc_symbol, function(x) x[which.max(x$mean_val)[1], ]))
rownames(best) <- best$hgnc_symbol
fpkm <- as.matrix(best[, samp_cols, drop = FALSE]); mode(fpkm) <- "numeric"
fpkm_log2 <- log2(fpkm + 0.5)
sg_158767 <- ifelse(grepl("^H", samp_cols), "H", "L"); names(sg_158767) <- samp_cols

for (g in c("H", "L")) {
  f <- file.path(diff158767_dir, paste0("GSE158767_limma_test_result.", g, "_vs_Others.txt"))
  vol_out <- file.path(diff158767_dir, "plots", paste0("GSE158767_volcano_", g, ".png"))
  hm_out  <- file.path(diff158767_dir, "plots", paste0("GSE158767_heatmap_", g, ".png"))
  if (file.exists(f)) {
    draw_volcano(f, paste0(g, " vs Others"), vol_out)
    sg <- ifelse(sg_158767 == g, g, "Others")
    draw_heatmap(fpkm_log2, sg, f, paste0("Heatmap: ", g, " vs Others (Top 50)"), hm_out)
  }
}

# ========================================================
# GSE99374 plots
# ========================================================
cat("\n=== GSE99374 ===\n")
diff99374_dir <- file.path(gse148036_dir, "差异化分析/GSE99374")
data99374_dir <- file.path(gse148036_dir, "GSE99374")

expr99374 <- read.table(file.path(data99374_dir, "G.txt"), header = TRUE, sep = "\t",
                         check.names = FALSE, quote = "", comment.char = "", row.names = 1)
cli99374 <- read.csv(file.path(data99374_dir, "GSE99374_cli.csv"), check.names = FALSE, row.names = 1)
subt99374 <- data.frame(Sample = cli99374[["title"]],
                        Group = cli99374[["disease group:ch1"]],
                        stringsAsFactors = FALSE)
subt99374 <- subt99374[subt99374$Sample %in% colnames(expr99374), ]
expr99374 <- as.matrix(expr99374[, subt99374$Sample, drop = FALSE]); mode(expr99374) <- "numeric"
expr99374_log2 <- log2(expr99374 + 1)

for (g in unique(subt99374$Group)) {
  sname <- safe_name(g)
  f <- file.path(diff99374_dir, paste0("GSE99374_limma_test_result.", sname, "_vs_Others.txt"))
  vol_out <- file.path(diff99374_dir, "plots", paste0("GSE99374_volcano_", sname, ".png"))
  hm_out  <- file.path(diff99374_dir, "plots", paste0("GSE99374_heatmap_", sname, ".png"))
  if (file.exists(f)) {
    draw_volcano(f, paste0(g, " vs Others"), vol_out)
    sg <- ifelse(subt99374$Group == g, g, "Others"); names(sg) <- subt99374$Sample
    draw_heatmap(expr99374_log2, sg, f, paste0("Heatmap: ", g, " vs Others (Top 50)"), hm_out)
  }
}

# ========================================================
# Immune infiltration heatmaps
# ========================================================
cat("\n=== Immune Infiltration Heatmaps ===\n")
immune_dir <- file.path(gse148036_dir, "免疫浸润")
for (dset in c("GSE148036", "GSE158767", "GSE99374")) {
  dir.create(file.path(immune_dir, dset, "plots"), showWarnings = FALSE, recursive = TRUE)
  tme_file <- file.path(immune_dir, dset, "tme_combine.csv")
  if (!file.exists(tme_file)) next

  dat <- read.csv(tme_file, check.names = FALSE)
  rownames(dat) <- dat$ID
  mat <- t(as.matrix(dat[, -1, drop = FALSE]))
  mode(mat) <- "numeric"

  pdf_path <- file.path(immune_dir, dset, "plots", "immune_heatmap.pdf")
  png_path <- file.path(immune_dir, dset, "plots", "immune_heatmap.png")

  pheatmap(mat, filename = png_path, width = 14, height = 10,
           color = colorRampPalette(c("#2b6cb0", "white", "#d94841"))(50),
           scale = "row", clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation", clustering_method = "ward.D2",
           fontsize = 8, fontsize_row = 7, show_colnames = TRUE,
           main = paste0(dset, " Immune Infiltration"))

  pheatmap(mat, filename = pdf_path, width = 14, height = 10,
           color = colorRampPalette(c("#2b6cb0", "white", "#d94841"))(50),
           scale = "row", clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation", clustering_method = "ward.D2",
           fontsize = 8, fontsize_row = 7, show_colnames = TRUE,
           main = paste0(dset, " Immune Infiltration"))

  cat(sprintf("  %s: immune_heatmap.png/pdf\n", dset))
}

cat("\nAll plots regenerated.\n")

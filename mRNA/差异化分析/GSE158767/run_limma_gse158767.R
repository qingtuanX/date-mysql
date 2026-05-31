suppressPackageStartupMessages({
  library(limma)
})

options(stringsAsFactors = FALSE)

data_dir <- normalizePath("../../GSE158767")
out_dir  <- "."
prefix   <- "GSE158767"
overwrite <- TRUE

raw_df <- read.csv(file.path(data_dir, "merged_data_dplyr.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)

sample_cols <- grep("^[HL][0-9]+$", colnames(raw_df), value = TRUE)
if (length(sample_cols) < 2) stop("No H/L sample columns found.")

if (!("hgnc_symbol" %in% colnames(raw_df))) stop("Missing hgnc_symbol column.")

raw_df$hgnc_symbol[is.na(raw_df$hgnc_symbol) | raw_df$hgnc_symbol == ""] <- NA
raw_df <- raw_df[!is.na(raw_df$hgnc_symbol), , drop = FALSE]

raw_df$mean_fpkm <- rowMeans(raw_df[, sample_cols, drop = FALSE], na.rm = TRUE)

best_per_gene <- do.call(rbind, by(raw_df, raw_df$hgnc_symbol, function(sub) {
  sub[which.max(sub$mean_fpkm)[1], , drop = FALSE]
}))
rownames(best_per_gene) <- best_per_gene$hgnc_symbol

expr_mat <- as.matrix(best_per_gene[, sample_cols, drop = FALSE])
mode(expr_mat) <- "numeric"

keep <- rowMeans(expr_mat) >= 1
expr_mat <- expr_mat[keep, , drop = FALSE]
cat(sprintf("Genes after mean FPKM >= 1 filter: %d\n", nrow(expr_mat)))

gset <- log2(expr_mat + 0.5)

subt <- data.frame(
  Sample = sample_cols,
  Group  = ifelse(grepl("^H", sample_cols), "H", "L"),
  stringsAsFactors = FALSE,
  row.names = sample_cols
)

run_one_vs_others <- function(log_expr, groups, target_group, output_file) {
  grp <- factor(ifelse(groups == target_group, target_group, "Others"),
                levels = c(target_group, "Others"))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)

  fit <- lmFit(log_expr, design)
  contrast <- makeContrasts(contrasts = paste0(target_group, " - Others"),
                            levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2, trend = TRUE)

  tt <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
  tt <- tt[, c("logFC", "t", "B", "P.Value", "adj.P.Val"), drop = FALSE]
  colnames(tt) <- c("log2FC", "t", "B", "PValue", "FDR")
  tt <- tt[order(tt$FDR), , drop = FALSE]

  write.table(tt, file = output_file, row.names = TRUE, col.names = NA,
              sep = "\t", quote = FALSE)
}

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

groups <- subt$Group
names(groups) <- subt$Sample
all_groups <- sort(unique(groups))

gset <- gset[, subt$Sample, drop = FALSE]

cat("Samples:", ncol(gset), "\n")
cat("Groups: ", paste(all_groups, collapse = " | "), "\n")
cat("Genes:  ", nrow(gset), "\n")

for (g in all_groups) {
  out_file <- file.path(out_dir,
    paste0(prefix, "_limma_test_result.", safe_name(g), "_vs_Others.txt"))
  if (file.exists(out_file) && !overwrite) {
    cat("Skipping existing:", out_file, "\n")
    next
  }
  run_one_vs_others(log_expr = gset,
                    groups = groups[colnames(gset)],
                    target_group = g,
                    output_file = out_file)
  cat("Output:", out_file, "\n")
}
cat("Analysis complete.\n")

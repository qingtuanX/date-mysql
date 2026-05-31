suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

options(stringsAsFactors = FALSE)

data_dir <- normalizePath("../../GSE99374")
out_dir  <- "."
prefix   <- "GSE99374"
overwrite <- TRUE

expr_file <- file.path(data_dir, "G.txt")
cli_file  <- file.path(data_dir, "GSE99374_cli.csv")

expr <- read.table(expr_file, header = TRUE, sep = "\t",
                   check.names = FALSE, quote = "", comment.char = "",
                   row.names = 1)
cli <- read.csv(cli_file, check.names = FALSE, stringsAsFactors = FALSE,
                row.names = 1)

if (!("title" %in% colnames(cli)) || !("disease group:ch1" %in% colnames(cli))) {
  stop("Clinical data missing required columns (title / disease group:ch1).")
}

subt <- data.frame(
  Sample = cli[["title"]],
  Group  = cli[["disease group:ch1"]],
  stringsAsFactors = FALSE
)
subt <- subt[!is.na(subt$Sample) & !is.na(subt$Group), , drop = FALSE]
subt <- subt[subt$Sample %in% colnames(expr), , drop = FALSE]

expr <- expr[, subt$Sample, drop = FALSE]
expr <- as.matrix(expr)
mode(expr) <- "numeric"

groups <- subt$Group
names(groups) <- subt$Sample

keep <- rowSums(cpm(expr) > 1) >= 2
expr <- expr[keep, , drop = FALSE]
cat(sprintf("Genes after CPM filter: %d / %d\n", nrow(expr), length(keep)))

run_one_vs_others <- function(count_mat, groups, target_group, output_file) {
  grp <- factor(ifelse(groups == target_group, "Target", "Others"),
                levels = c("Target", "Others"))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- c("Target", "Others")

  dge <- DGEList(counts = count_mat)
  dge <- calcNormFactors(dge)
  v   <- voom(dge, design)

  fit <- lmFit(v, design)
  contrast <- makeContrasts(contrasts = "Target - Others",
                            levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)

  tt <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
  tt <- tt[, c("logFC", "t", "B", "P.Value", "adj.P.Val"), drop = FALSE]
  colnames(tt) <- c("log2FC", "t", "B", "PValue", "FDR")
  tt <- tt[order(tt$FDR), , drop = FALSE]

  write.table(tt, file = output_file, row.names = TRUE, col.names = NA,
              sep = "\t", quote = FALSE)
}

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

all_groups <- sort(unique(groups))
cat("Samples:", ncol(expr), "\n")
cat("Groups: ", paste(all_groups, collapse = " | "), "\n")
cat("Genes:  ", nrow(expr), "\n")

for (g in all_groups) {
  out_file <- file.path(out_dir,
    paste0(prefix, "_limma_test_result.", safe_name(g), "_vs_Others.txt"))
  if (file.exists(out_file) && !overwrite) {
    cat("Skipping existing:", out_file, "\n")
    next
  }
  run_one_vs_others(count_mat = expr,
                    groups = groups[colnames(expr)],
                    target_group = g,
                    output_file = out_file)
  cat("Output:", out_file, "\n")
}
cat("Analysis complete.\n")

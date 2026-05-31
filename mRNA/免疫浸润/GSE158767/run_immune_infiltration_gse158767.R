# ================================================================
# GSE158767 免疫浸润分析脚本
# 输入: ENST transcript ID -> SYMBOL -> IOBR
# ================================================================

rm(list = ls())
setwd('d:/CODE/date-mysql/mRNA/GSE158767/')
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

cat("========================================\n")
cat("GSE158767 免疫浸润分析 - 开始\n")
cat("========================================\n")

# 依赖包
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
if (!require("readr", quietly = TRUE)) install.packages("readr")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!require("GSVA", quietly = TRUE)) BiocManager::install("GSVA")
if (!require("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!require("sva", quietly = TRUE)) BiocManager::install("sva")
lock_path <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.5/00LOCK-org.Hs.eg.db")
if (dir.exists(lock_path)) unlink(lock_path, recursive = TRUE, force = TRUE)
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

if (!require("IOBR", quietly = TRUE)) {
  remotes::install_github("IOBR/IOBR", upgrade = "never")
}

library(IOBR)
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

if (!exists("deconvo_tme", mode = "function")) {
  stop("IOBR已加载，但未找到deconvo_tme函数")
}

cat("\n[步骤 1] 读取表达矩阵...\n")
exp_raw <- read.csv("GSE158767_mRNA.fpkm.csv", check.names = FALSE)
gene_col <- as.character(exp_raw[, 1])
exp_raw <- exp_raw[, -1, drop = FALSE]

cat("原始维度: ", nrow(exp_raw), " x ", ncol(exp_raw), "\n", sep = "")

# ENSTxxxx.xx -> ENSTxxxx
ids <- gsub("\\..*$", "", gene_col)

cat("\n[步骤 2] 转换ID: ENST -> SYMBOL...\n")
symbols <- mapIds(
  org.Hs.eg.db,
  keys = ids,
  keytype = "ENSEMBLTRANS",
  column = "SYMBOL",
  multiVals = "first"
)

mapped <- !is.na(symbols) & symbols != ""
cat("可映射条目: ", sum(mapped), " / ", length(symbols), "\n", sep = "")

exp_mapped <- exp_raw[mapped, , drop = FALSE]
exp_mapped$SYMBOL <- symbols[mapped]

# 按SYMBOL聚合（取均值）
exp_symbol <- exp_mapped %>%
  as.data.frame() %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

exp_symbol$SYMBOL <- toupper(exp_symbol$SYMBOL)
rownames(exp_symbol) <- exp_symbol$SYMBOL
exp_symbol <- exp_symbol[, -1, drop = FALSE]

eset <- as.matrix(exp_symbol)
mode(eset) <- "numeric"

cat("SYMBOL矩阵维度: ", nrow(eset), " x ", ncol(eset), "\n", sep = "")

results_list <- list()
array <- FALSE

safe_run <- function(tag, expr) {
  cat("  - 运行 ", tag, "...\n", sep = "")
  tryCatch(
    expr,
    error = function(e) {
      msg <- conditionMessage(e)
      if (is.null(msg) || msg == "") msg <- "未知错误"
      cat("    ✗ ", tag, " 失败: ", msg, "\n", sep = "")
      NULL
    }
  )
}

cat("\n[步骤 3] 免疫浸润分析...\n")
cibersort <- safe_run("CIBERSORT", deconvo_tme(eset = eset, method = "cibersort", arrays = array, perm = 100))
if (!is.null(cibersort)) results_list$cibersort <- cibersort

epic <- safe_run("EPIC", deconvo_tme(eset = eset, method = "epic", arrays = array))
if (!is.null(epic)) results_list$epic <- epic

mcp <- safe_run("MCPcounter", deconvo_tme(eset = eset, method = "mcpcounter"))
if (!is.null(mcp)) results_list$mcp <- mcp

xcell <- safe_run("xCell", deconvo_tme(eset = eset, method = "xcell", arrays = array))
if (!is.null(xcell)) results_list$xcell <- xcell

estimate <- safe_run("ESTIMATE", deconvo_tme(eset = eset, method = "estimate"))
if (!is.null(estimate)) results_list$estimate <- estimate

timer <- safe_run("TIMER", deconvo_tme(eset = eset, method = "timer", group_list = rep("dlbc", ncol(eset))))
if (!is.null(timer)) results_list$timer <- timer

quantiseq <- safe_run("quanTIseq", deconvo_tme(eset = eset, method = "quantiseq", tumor = TRUE, arrays = array, scale_mrna = TRUE))
if (!is.null(quantiseq)) results_list$quantiseq <- quantiseq

ips <- safe_run("IPS", deconvo_tme(eset = eset, method = "ips", plot = FALSE))
if (!is.null(ips)) results_list$ips <- ips

cat("\n[步骤 4] 保存结果...\n")
if (length(results_list) == 0) {
  stop("没有任何算法成功，无法生成结果。")
}

tme_combine <- results_list[[1]]
if (length(results_list) > 1) {
  for (i in 2:length(results_list)) {
    tme_combine <- inner_join(tme_combine, results_list[[i]], by = "ID")
  }
}

if (exists("cibersort")) save(cibersort, file = "cibersort-TME-Cell-fration.RData")
if (exists("epic")) save(epic, file = "epic-TME-Cell-fration.RData")
if (exists("mcp")) save(mcp, file = "mcp-TME-Cell-fration.RData")
if (exists("xcell")) save(xcell, file = "xcell-TME-Cell-fration.RData")
if (exists("estimate")) save(estimate, file = "estimate-TME-Cell-fration.RData")
if (exists("timer")) save(timer, file = "timer-TME-Cell-fration.RData")
if (exists("quantiseq")) save(quantiseq, file = "quantiseq-TME-Cell-fration.RData")
if (exists("ips")) save(ips, file = "ips-TME-Cell-fration.RData")
save(tme_combine, file = "tme_combine.RData")

if (exists("cibersort")) write.csv(cibersort, file = "cibersort.csv", row.names = FALSE)
if (exists("epic")) write.csv(epic, file = "epic.csv", row.names = FALSE)
if (exists("mcp")) write.csv(mcp, file = "mcp.csv", row.names = FALSE)
if (exists("xcell")) write.csv(xcell, file = "xcell.csv", row.names = FALSE)
if (exists("estimate")) write.csv(estimate, file = "estimate.csv", row.names = FALSE)
if (exists("timer")) write.csv(timer, file = "timer.csv", row.names = FALSE)
if (exists("quantiseq")) write.csv(quantiseq, file = "quantiseq.csv", row.names = FALSE)
if (exists("ips")) write.csv(ips, file = "ips.csv", row.names = FALSE)
write.csv(tme_combine, file = "tme_combine.csv", row.names = FALSE)

cat("✓ 分析完成\n")
cat("输出目录: ", getwd(), "\n", sep = "")

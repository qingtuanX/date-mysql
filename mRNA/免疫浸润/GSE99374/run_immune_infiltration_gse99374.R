# ================================================================
# GSE99374 免疫浸润分析脚本
# 基于IOBR包，进行多种免疫细胞浸润推断
# ================================================================

rm(list = ls())
setwd('d:/CODE/date-mysql/mRNA/GSE99374/')

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

cat("========================================\n")
cat("GSE99374 免疫浸润分析 - 开始\n")
cat("========================================\n")

# ================================================================
# 步骤0: 安装必要的包
# ================================================================
cat("\n[步骤 0] 安装依赖包...\n")

# 安装BiocManager (如果需要)
if (!require("BiocManager", quietly = TRUE)) {
  cat("安装 BiocManager...\n")
  install.packages("BiocManager")
}

# 安装devtools
if (!require("devtools", quietly = TRUE)) {
  cat("安装 devtools...\n")
  install.packages("devtools")
}

# 安装remotes（devtools::install_github底层依赖）
if (!require("remotes", quietly = TRUE)) {
  cat("安装 remotes...\n")
  install.packages("remotes")
}

# 安装GSVA
if (!require("GSVA", quietly = TRUE)) {
  cat("安装 GSVA...\n")
  BiocManager::install("GSVA")
}

# 安装readxl
if (!require("readxl", quietly = TRUE)) {
  cat("安装 readxl...\n")
  install.packages("readxl")
}

cat("✓ 依赖包安装完成\n")

# ================================================================
# 步骤1: 安装和加载IOBR包
# ================================================================
cat("\n[步骤 1] 安装IOBR包...\n")

if (!require("IOBR", quietly = TRUE)) {
  cat("从GitHub安装IOBR...\n")
  tryCatch({
    remotes::install_github("IOBR/IOBR", upgrade = "never")
    library(IOBR)
    cat("✓ IOBR安装成功\n")
  }, error = function(e) {
    cat("✗ IOBR安装失败: ", e$message, "\n", sep = "")
  })
} else {
  cat("✓ IOBR已安装\n")
}

# 强校验：未成功加载IOBR时直接停止，避免后续deconvo_tme连续失败
if (!require("IOBR", quietly = TRUE)) {
  stop("IOBR未成功安装/加载。请先安装IOBR后重试。")
}

# 强校验：函数可用性
if (!exists("deconvo_tme", mode = "function")) {
  stop("IOBR已加载，但未找到deconvo_tme函数。请重装IOBR并重试。")
}

# 加载其他必要的库
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(readxl, quietly = TRUE)

# ================================================================
# 步骤2: 载入基因表达数据
# ================================================================
cat("\n[步骤 2] 加载基因表达数据...\n")

# 优先使用G.txt（基因符号矩阵，适合IOBR算法），否则回退到GSE99374_exp.csv
if (file.exists("G.txt")) {
  eset_raw <- read.delim("G.txt", check.names = FALSE)
  eset <- as.data.frame(eset_raw)
  rownames(eset) <- eset[, 1]
  eset <- eset[, -1, drop = FALSE]
  cat("✓ 使用 G.txt 作为输入矩阵\n")
} else {
  eset_raw <- read_excel("GSE99374_exp.csv")
  eset <- as.data.frame(eset_raw)
  rownames(eset) <- eset[, 1]
  eset <- eset[, -1, drop = FALSE]
  cat("✓ 使用 GSE99374_exp.csv 作为输入矩阵\n")
}

# 去重与缺失处理
eset <- eset[!duplicated(rownames(eset)), , drop = FALSE]
eset <- eset[rownames(eset) != "" & !is.na(rownames(eset)), , drop = FALSE]

# 转换为数值矩阵
eset <- as.matrix(eset)
mode(eset) <- "numeric"

cat("✓ 表达数据加载成功\n")
cat("  维度: ", nrow(eset), " 基因 x ", ncol(eset), " 样本\n")
cat("  样本: ", paste(colnames(eset), collapse = ", "), "\n")

# ================================================================
# 步骤3: 使用多种算法进行免疫细胞浸润推断
# ================================================================
cat("\n[步骤 3] 进行免疫细胞浸润分析...\n")

array <- FALSE  # FALSE表示RNA-seq数据
results_list <- list()

# CIBERSORT方法
cat("  - 运行 CIBERSORT...\n")
tryCatch({
  cibersort <- deconvo_tme(eset = eset, method = "cibersort", arrays = array, perm = 100)
  results_list$cibersort <- cibersort
  cat("    ✓ CIBERSORT 完成 (", nrow(cibersort), " x ", ncol(cibersort), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ CIBERSORT 失败:", e$message, "\n")
})

# EPIC方法
cat("  - 运行 EPIC...\n")
tryCatch({
  epic <- deconvo_tme(eset = eset, method = "epic", arrays = array)
  results_list$epic <- epic
  cat("    ✓ EPIC 完成 (", nrow(epic), " x ", ncol(epic), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ EPIC 失败:", e$message, "\n")
})

# MCP计数方法
cat("  - 运行 MCPcounter...\n")
tryCatch({
  mcp <- deconvo_tme(eset = eset, method = "mcpcounter")
  results_list$mcp <- mcp
  cat("    ✓ MCPcounter 完成 (", nrow(mcp), " x ", ncol(mcp), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ MCPcounter 失败:", e$message, "\n")
})

# xCell方法
cat("  - 运行 xCell...\n")
tryCatch({
  xcell <- deconvo_tme(eset = eset, method = "xcell", arrays = array)
  results_list$xcell <- xcell
  cat("    ✓ xCell 完成 (", nrow(xcell), " x ", ncol(xcell), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ xCell 失败:", e$message, "\n")
})

# ESTIMATE方法
cat("  - 运行 ESTIMATE...\n")
tryCatch({
  estimate <- deconvo_tme(eset = eset, method = "estimate")
  results_list$estimate <- estimate
  cat("    ✓ ESTIMATE 完成 (", nrow(estimate), " x ", ncol(estimate), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ ESTIMATE 失败:", e$message, "\n")
})

# TIMER方法
cat("  - 运行 TIMER...\n")
tryCatch({
  # TIMER要求提供已知肿瘤类型；此处使用dlbc作为占位类型以满足接口要求
  timer <- deconvo_tme(eset = eset, method = "timer", group_list = rep("dlbc", dim(eset)[2]))
  results_list$timer <- timer
  cat("    ✓ TIMER 完成 (", nrow(timer), " x ", ncol(timer), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ TIMER 失败:", e$message, "\n")
})

# quanTIseq方法
cat("  - 运行 quanTIseq...\n")
tryCatch({
  quantiseq <- deconvo_tme(eset = eset, method = "quantiseq", tumor = TRUE, arrays = array, scale_mrna = TRUE)
  results_list$quantiseq <- quantiseq
  cat("    ✓ quanTIseq 完成 (", nrow(quantiseq), " x ", ncol(quantiseq), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ quanTIseq 失败:", e$message, "\n")
})

# IPS方法
cat("  - 运行 IPS...\n")
tryCatch({
  ips <- deconvo_tme(eset = eset, method = "ips", plot = FALSE)
  results_list$ips <- ips
  cat("    ✓ IPS 完成 (", nrow(ips), " x ", ncol(ips), ")\n", sep="")
}, error = function(e) {
  cat("    ✗ IPS 失败:", e$message, "\n")
})

# ================================================================
# 步骤4: 合并所有结果
# ================================================================
cat("\n[步骤 4] 合并所有算法结果...\n")

if (length(results_list) > 0) {
  tme_combine <- results_list[[1]]
  
  if (length(results_list) > 1) {
    for (i in 2:length(results_list)) {
      tme_combine <- inner_join(tme_combine, results_list[[i]], by = "ID")
    }
  }
  
  cat("✓ 合并完成\n")
  cat("  合并后维度: ", nrow(tme_combine), " x ", ncol(tme_combine), "\n", sep="")
  
  # ================================================================
  # 步骤5: 保存RData格式
  # ================================================================
  cat("\n[步骤 5] 保存RData文件...\n")
  
  if (exists("cibersort")) save(cibersort, file = "cibersort-TME-Cell-fration.RData")
  if (exists("epic")) save(epic, file = "epic-TME-Cell-fration.RData")
  if (exists("mcp")) save(mcp, file = "mcp-TME-Cell-fration.RData")
  if (exists("xcell")) save(xcell, file = "xcell-TME-Cell-fration.RData")
  if (exists("estimate")) save(estimate, file = "estimate-TME-Cell-fration.RData")
  if (exists("timer")) save(timer, file = "timer-TME-Cell-fration.RData")
  if (exists("quantiseq")) save(quantiseq, file = "quantiseq-TME-Cell-fration.RData")
  if (exists("ips")) save(ips, file = "ips-TME-Cell-fration.RData")
  save(tme_combine, file = "tme_combine.RData")
  
  cat("✓ RData文件保存完成\n")
  
  # ================================================================
  # 步骤6: 导出CSV格式
  # ================================================================
  cat("\n[步骤 6] 保存CSV文件...\n")
  
  if (exists("cibersort")) write.csv(cibersort, file = "cibersort.csv", row.names = FALSE)
  if (exists("mcp")) write.csv(mcp, file = "mcp.csv", row.names = FALSE)
  if (exists("xcell")) write.csv(xcell, file = "xcell.csv", row.names = FALSE)
  if (exists("estimate")) write.csv(estimate, file = "estimate.csv", row.names = FALSE)
  if (exists("epic")) write.csv(epic, file = "epic.csv", row.names = FALSE)
  if (exists("quantiseq")) write.csv(quantiseq, file = "quantiseq.csv", row.names = FALSE)
  if (exists("ips")) write.csv(ips, file = "ips.csv", row.names = FALSE)
  if (exists("timer")) write.csv(timer, file = "timer.csv", row.names = FALSE)
  write.csv(tme_combine, file = "tme_combine.csv", row.names = FALSE)
  
  cat("✓ CSV文件保存完成\n")
  
} else {
  cat("✗ 未能成功运行任何分析方法\n")
}

cat("\n========================================\n")
cat("应吧 分析完成！\n")
cat("========================================\n")
cat("\n输出文件位置: ", getwd(), "\n", sep="")

# 列出生成的文件
output_files <- list.files(pattern = "\\.(csv|RData)$")
if (length(output_files) > 0) {
  cat("\n生成的文件:\n")
  for (f in output_files) {
    cat("  - ", f, "\n", sep="")
  }
}

cat("\n")

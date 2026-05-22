# ================================================================
# GSE99374 免疫浸润分析 - 热力图可视化脚本
# 使用tme_combine.csv结果制作热力图
# ================================================================

rm(list = ls())
setwd('d:/CODE/date-mysql/mRNA/GSE99374/')

# 设置CRAN镜像，避免Windows下install.packages无镜像报错
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

cat("========================================\n")
cat("免疫浸润热力图生成 - 开始\n")
cat("========================================\n")

# ================================================================
# 步骤 1: 安装和加载必要的包
# ================================================================
cat("\n[步骤 1] 加载包...\n")

required_packages <- c("pheatmap", "RColorBrewer")

for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("✓ 包加载成功\n")

# ================================================================
# 步骤 2: 加载数据
# ================================================================
cat("\n[步骤 2] 加载免疫浸润数据...\n")

# 检查文件
if (!file.exists("tme_combine.csv")) {
  cat("✗ 错误：找不到tme_combine.csv文件\n")
  cat("  请先运行: run_immune_infiltration_gse99374_v2.R\n")
  stop("数据文件不存在")
}

# 加载数据
tme_data <- read.csv("tme_combine.csv")
cat("✓ 数据加载成功\n")
cat("  维度:", nrow(tme_data), "样本 x", ncol(tme_data), "列\n")

# 提取样本ID和免疫细胞数据
sample_ids <- tme_data$ID
immune_data <- tme_data[, -1]
rownames(immune_data) <- sample_ids

cat("  样本名称:", paste(sample_ids, collapse = ", "), "\n")

# ================================================================
# 步骤 3: 创建热力图
# ================================================================
cat("\n[步骤 3] 生成热力图...\n")

# 将数据转置（样本为列，细胞类型为行）
heatmap_data <- t(immune_data)

# 创建PDF文件
pdf("免疫浸润热力图.pdf", width = 12, height = 10)

cat("  生成: 免疫浸润热力图.pdf\n")

# 创建热力图
pheatmap(heatmap_data,
         main = "GSE99374 - 免疫细胞浸润热力图",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         scale = "row",  # 按行标准化
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         fontsize = 10,
         fontsize_row = 8,
         show_colnames = TRUE,
         show_rownames = TRUE,
         legend = TRUE)

dev.off()
cat("  ✓ 已保存\n")

# ================================================================
# 步骤 4: 生成总结
# ================================================================
cat("\n[步骤 4] 生成分析总结...\n")

summary_text <- paste(
  "GSE99374 免疫浸润热力图分析\n",
  "=====================================\n\n",
  
  "样本信息:\n",
  "  - 总样本数: ", length(sample_ids), "\n",
  "  - 样本名称: ", paste(sample_ids, collapse = ", "), "\n\n",
  
  "免疫细胞类型:\n",
  "  - 总细胞类型数: ", ncol(immune_data), "\n",
  "  - 细胞类型列表:\n",
  paste("    * ", colnames(immune_data)[1:min(15, ncol(immune_data))], collapse = "\n    * "), "\n",
  ifelse(ncol(immune_data) > 15, paste("\n    ... 及其他", ncol(immune_data) - 15, "种"), ""), "\n\n",
  
  "热力图说明:\n",
  "  - 红色: 较高的免疫细胞浸润值\n",
  "  - 白色: 中等的免疫细胞浸润值\n",
  "  - 蓝色: 较低的免疫细胞浸润值\n",
  "  - 树状图: 显示相似样本和细胞类型聚集\n",
  "  - 行标准化: 使不同细胞类型的浸润水平可比较\n\n",
  
  "如何解读热力图:\n",
  "  1. 观察颜色模式 - 红色聚集区域表示该样本中该细胞类型浸润高\n",
  "  2. 查看树状图 - 距离近的样本免疫特征相似\n",
  "  3. 比较样本分组:\n",
  "     - TB组(TB1-TB3): 结核病患者\n",
  "     - CD组(CD1-CD3): 对照组\n",
  "     - N组(N-2-N-4): 健康人群\n",
  "  4. 识别差异模式 - 不同分组是否显示不同的免疫浸润模式\n\n",
  
  "生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  sep = ""
)

writeLines(summary_text, "热力图分析说明.txt")
cat(summary_text)

# ================================================================
# 完成
# ================================================================
cat("\n========================================\n")
cat("✓ 热力图生成完成！\n")
cat("========================================\n")
cat("\n生成的文件:\n")
cat("  1. 免疫浸润热力图.pdf\n")
cat("  2. 热力图分析说明.txt\n")
cat("\n文件位置: d:/CODE/date-mysql/mRNA/GSE99374/\n\n")

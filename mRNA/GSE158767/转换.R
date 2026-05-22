rm(list=ls())
setwd("D:\\文章\\文章投稿\\结核病宿主数据库构建\\mRNA\\GSE158767")
library(biomaRt)

# 连接到ENSEMBL数据库
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
exp<-read.csv(file="GSE158767_mRNA.fpkm.csv")
# 准备你的转录本ID
transcript_ids <- exp$Transcript_id

# 获取映射信息
results <- getBM(
  attributes = c("ensembl_transcript_id_version", "hgnc_symbol", 
                 "ensembl_gene_id", "transcript_biotype"),
  filters = "ensembl_transcript_id_version",
  values = transcript_ids,
  mart = ensembl
)

# 查看结果
print(results)

# 如果想要直接对应关系
mapping <- data.frame(
  transcript_id = transcript_ids,
  gene_symbol = sapply(transcript_ids, function(x) {
    sym <- results$hgnc_symbol[results$ensembl_transcript_id_version == x]
    if (length(sym) > 0) sym[1] else NA
  })
)
print(mapping)



library(dplyr)
colnames(results)[colnames(results) == "ensembl_transcript_id_version"] <- "Transcript_id"
# 1. 内连接
inner_join_df <- inner_join(exp, results, by = "Transcript_id")

# 如果有相同列名（除了Transcript_id），dplyr会自动添加后缀 .x 和 .y
# 可以手动重命名
merged_df <- inner_join(exp, results, by = "Transcript_id", suffix = c("_file1", "_file2"))

# 保存
write.csv(merged_df, "merged_data_dplyr.csv", row.names = FALSE)





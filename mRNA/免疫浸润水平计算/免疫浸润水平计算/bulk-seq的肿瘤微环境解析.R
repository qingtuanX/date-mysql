
rm (list = ls ())
# begin
if (!require("devtools")) 
  install.packages("devtools", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(download.file.method = "wininet")

# 加载 devtools
library(devtools)
if (!requireNamespace("IOBR", quietly = TRUE)) 
 devtools::install_github("IOBR/IOBR")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")

setwd("C:\\Users\\19924\\Desktop\\r_iobr")
list.files()
library(IOBR)
library(dplyr)
library (tidyr)
library(survminer)

eset=read.table("gene_expression.txt",row.names=1,header=T)
eset[1:5, 1:5]



epic     <-deconvo_tme(eset = eset,method = "epic",arrays = array)
mcp      <-deconvo_tme(eset = eset,method = "mcpcounter")
xcell    <-deconvo_tme(eset = eset,method = "xcell",arrays = array)
estimate <-deconvo_tme(eset = eset,method = "estimate")

#timer_available_cancers
timer    <-deconvo_tme(eset = eset,method = "timer",group_list = rep("dlbc",dim(eset)[2]))

quantiseq<-deconvo_tme(eset = eset,method = "quantiseq", tumor = TRUE, arrays = array, scale_mrna = TRUE)
res_quantiseq<-cell_bar_plot(input = quantiseq, id = "ID", title = "quanTIseq Cell Fraction")

ips      <-deconvo_tme(eset = eset,method = "ips",plot= FALSE)

#?????????߽??󷽷??ϲ???һ??
tme_combine<-cibersort %>% 
  inner_join(.,mcp,by       = "ID") %>% 
  inner_join(.,xcell,by     = "ID") %>%
  inner_join(.,epic,by      = "ID") %>% 
  inner_join(.,estimate,by  = "ID") %>% 
  inner_join(.,timer,by     = "ID") %>% 
  inner_join(.,quantiseq,by = "ID") %>% 
  inner_join(.,ips,by       = "ID")
dim(tme_combine)

save(cibersort,file = "cibersort-TME-Cell-fration.RData")
save(epic,file = "epic-TME-Cell-fration.RData")
save(mcp,file = "mcp-TME-Cell-fration.RData")
save(xcell,file = "xcell-TME-Cell-fration.RData")
save(estimate,file = "estimate-TME-Cell-fration.RData")
save(timer,file = "timer-TME-Cell-fration.RData")
save(quantiseq,file = "quantiseq-TME-Cell-fration.RData")
save(ips,file = "ips-TME-Cell-fration.RData")
save(tme_combine,file = "tme_combine.RData")

#??cibersort??xcell??estimate??epicΪ???????????洢??csv?У??????㷨????
write.csv(cibersort, file= "cibe112rsort.csv", row.names = F)
write.csv(xcell, file= "xcell.csv", row.names = F)
write.csv(estimate, file= "esti1mate.csv", row.names = F)
write.csv(epic, file= "epic.csv", row.names = F)
write.csv(quantiseq, file= "estimate.csv", row.names = F)
write.csv(ips, file= "ips.csv", row.names = F)
write.csv(timer, file= "timer.csv", row.names = F)
write.csv(tme_combine, file= "tme_combine.csv", row.names = F)



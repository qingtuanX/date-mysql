rm(list=ls())#clear Global Environment
setwd('D:/桌面/SCI论文写作与绘图/R语言绘图/基础图形绘制/多组火山图')#设置工作路径

#加载R包
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales) # Scale Functions for Visualization

#加载数据（随机编写，无实际意义）
df <- read.table("data.txt", header = 1, check.names = F, sep = "\t")
df$group <- factor(df$group, levels = c("group1-group2", "group1-group3", "group2-group3"))
##与之前绘制单组火山图一致，先根据设定阈值确定所有OTU的显著性
df$group2<-as.factor(ifelse(df$p_value < 0.05 & abs(df$log2FC) >= 2, 
                           ifelse(df$log2FC>= 2 ,'up','down'),'NS'))
df$group2 <- factor(df$group2, levels = c("up", "down", "NS"))
##确定添加标签的数据
df$label<-ifelse(df$p_value<0.05&abs(df$log2FC)>=4,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$OTU), '')

##为了构建图形中所有数据点背景框，需要先确定每个组的最大值与最小值
df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(log2FC),min_log2FC = min(log2FC))

##绘制图中各组数据点中的灰色背景柱子
p <- ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)
p

##添加各组的数据点
p1 <- p+geom_jitter(data = df,
                   mapping = aes(x = group, y = log2FC, color = group2),
                   size= 3,width =0.4, alpha = 0.7)
p1

##将分组信息通过在X=0的位置添加方块进行展示
#有两种方式可以实现，1）和灰色背景柱子添加方式一致，采用geom_col函数添加
p2 <- p1+geom_col(data = df_bg,
                  mapping = aes(x= group, y = 0.3, fill = group))+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -0.3, fill = group))
p2
#2）通过geom_rect函数手动添加,分组太多不推荐使用
p2 <- p1+geom_rect(data = df_bg,
                   aes(xmin = 1-0.4, xmax = 1+0.4, ymin = -0.3, ymax = 0.3),
                   alpha = 0.5,
                   color = NA,
                   fill = "red",
                   show.legend = F)+
  geom_rect(data = df_bg,
            aes(xmin = 2-0.4, xmax = 2+0.4, ymin = -0.3, ymax = 0.3),
            alpha = 0.5,
            color = NA,
            fill = "green",
            show.legend = F)+
  geom_rect(data = df_bg,
            aes(xmin = 3-0.4, xmax = 3+0.4, ymin = -0.3, ymax = 0.3),
            alpha = 0.5,
            color = NA,
            fill = "yellow",
            show.legend = F)
  
p2

##在方块中添加分组的文字信息
p3 <- p2+geom_text(data=df_bg,
                   mapping = aes(x=group, y=0, label=group),
                   size = 4, color ="#dbebfa")
p3


####个性化绘图模板
ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  #添加各组的数据点
  geom_jitter(data = df,
              mapping = aes(x = group, y = log2FC, color = group2),
              size= 3,width = 0.4, alpha = 0.7)+
  # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
  geom_col(data = df_bg,
           mapping = aes(x= group, y = 0.4, fill = group),
           width = 0.8)+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -0.4, fill = group),
           width = 0.8)+
  # 在方块中添加分组的文字信息
  geom_text(data=df_bg,
            mapping = aes(x=group, y=0, label=group),
            size = 4, color ="#dbebfa",fontface = "bold")+
  #根据需要决定是添加辅助线
  # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
  # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
  #颜色设置
  scale_color_manual(values = c("#e42313", "#0061d5", "#8b8c8d"))+
  scale_fill_manual(values = c("#ed7902", "#ef5734", "#b5c327"))+
  #添加显著性标签
  geom_text_repel(data = df,
                  mapping = aes(x = group, y = log2FC, label = label),
                  max.overlaps = 10000,
                  size=3,
                  box.padding=unit(0.8,'lines'),
                  point.padding=unit(0.8, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)+
  # 主题设置
  theme_classic()+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.8))+
  labs(x = "group", y = "Log2FoldChange", fill= NULL, color = NULL)+
  #调整图例
  guides(color=guide_legend(override.aes = list(size=6,alpha=1)))

#背景色
color <- colorRampPalette(brewer.pal(11,"BrBG"))(30)
#添加背景
grid.raster(alpha(color, 0.2), 
            width = unit(1, "npc"), 
            height = unit(1,"npc"),
            interpolate = T)


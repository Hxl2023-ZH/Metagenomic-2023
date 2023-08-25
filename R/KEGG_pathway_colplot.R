#hxl

setwd("C:/Users/huxin/Desktop/Final_result_analysys/Bacteria_Emapper_annotations_03/KEGG/Pathway/lefse0.001/")

library(ggplot2)
library(ggpubr)

df <- read.delim("pathwaylefse_0.001.txt", check.names = F)
head(df)
#df <- df %>% mutate(Group = factor(Group, levels = c("Obesity", "Control")))
#df$Group <- as.factor(df$Group, levels=c("Obesity", "Control"))
df <- df %>% mutate(Features = factor(Features, levels = Features))
df$Group <- as.factor(df$Group)
head(df)
####使用ggplot绘制条形图
ggplot(df, aes(x=LDA, y=Features, fill=Group)) + #构建ggplot绘图几何对象
  geom_col(width = 0.9) +#绘制条形图并设置大小
  scale_fill_manual(values = c('#CD5B45','#00688B'))+#手动设置颜色
  theme_minimal()+#设置主题
  theme(panel.grid =element_blank())+#去除网格线
  #theme(panel.border = element_blank())+#去除边框
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        legend.position = c( 0.93, 0.93),
        #legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size =13))+
  geom_text(data=subset(df, LDA>0), #设置正值的注释文本的位置和大小
            aes(x=0, y= Features,label= paste0(Features, ' ')), 
            size = 4, hjust = 'outward' ) +
  geom_text(data=subset(df, LDA<0), #设置负值的注释文本的位置和大小
            aes(x=0, y= Features,label= paste0(' ', Features)), 
            size = 4, hjust = 'inward' ) +
  geom_vline(xintercept = c(-2, -3, 0, 2, 3, 4),lty="dashed")+#添加垂直线
  labs(x = "LDA SCORE (log10)", y = "")#设置X,Y轴标题文本

ggsave(filename = "all-7_kegg_pathway_lefse_001.png", width = 12, height = 5, units = "in", dpi = 300)
ggsave(filename = "all-7_kegg_pathway_lefse_001.pdf", width = 12, height = 7)


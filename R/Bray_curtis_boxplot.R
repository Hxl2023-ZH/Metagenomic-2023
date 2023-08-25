
################################################################################
###                             Beta多样性分析
################################################################################
# XinLiang Hu
# 2023-08-23

#设置工作目录
setwd("C:/Users/huxin/Desktop/Final_result_analysys/01_Bacteria_species_abundance/Belta-diversity/")

library(ggplot2)

#加载数据表
df <- read.delim('total-7_abundance_table_species_0.txt', row.names = 1,check.names=F)
group <- read.delim('total-7_metadata1.txt', stringsAsFactors = FALSE,check.names=F)

#计算Beta多样性群落相异指数
dis <- vegan::vegdist(df, method = 'bray')

dis <- as.matrix(dis)
#class(dis)
#[1] "matrix" "array"
write.table(dis, 'total7_bray-curtis-0.txt', sep = '\t', col.names = NA, quote = FALSE)

dis <- read.delim('total7_bray-curtis-0.txt', row.names = 1,check.names=F)
#class(dis)
#[1] "data.frame"

#根据分组获得组内距离矩阵
Control <- subset(group, Type == 'Control')$Individuals
head(Control)
dis_Control <- dis[Control,Control]
head(dis_Control)

Obesity <- subset(group, Type == 'Obesity')$Individuals
dis_Obesity <- dis[Obesity,Obesity]
head(dis_Obesity)

#将矩阵转化为向量，以便用于作图和统计
dis_Control <- as.vector(as.dist(dis_Control))
dis_Obesity <- as.vector(as.dist(dis_Obesity))


#构建作图数据集
dat <- data.frame(
  dis = c(dis_Control, dis_Obesity),
  group = factor(c(
    rep('Control', length(dis_Control)), 
    rep('Obesity', length(dis_Obesity))
  ), levels = c('Control', 'Obesity'))
)

head(dat)


#绘制 Bray-curtis 距离指数分布的箱线图
ggplot(dat, aes(group, dis)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_manual(values = c('#CD5B45','#00688B'))+
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        legend.position = 'none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15)) +
  geom_signif(comparisons = list(c("Control", "Obesity")), 
              test = wilcox.test, 
              map_signif_level = T, 
              color = "black", textsize = 6) +
  labs(x = "", y = 'Bray-Curtis Dissimilarity')


ggsave(filename = "all-7_beta-bray_box.png", width = 5, height = 7, units = "in", dpi = 300)
ggsave(filename = "all-7_beta-bray_box.pdf", width = 5, height = 7)


############## 以下部分可以不运行 ###################
#使用 Kruskal-Wallis Test，组间的整体差异检验
kruskal.test(dis~group, data = dat)

wilcox.test(dis_Control, dis_Obesity, alternative = 'two.sided')

#考虑到 Wilcoxon 秩和检验体现了中位数的差异，因此计算三组数据的中位数以评估 Beta 多样性的高低水平
median(dis_Control)
#[1] 0.7686048
median(dis_Obesity)
#[1] 0.7844391

#基于上述统计结果，判断好组间差异后，将差异分析结果添加到箱线图中
p +
  annotate('text', label = 'Kruskal-Wallis Test', x = 1.5, y = 0.95, size = 3) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.001), x = 1.5, y = 0.92, size = 3, parse = TRUE)

ggsave(filename = "all7_beta-bray_box.png", width = 7, height = 7, units = "in", dpi = 300)

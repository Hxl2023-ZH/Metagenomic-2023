
################################################################################
###                      Alpha多样性分析                                    ###
################################################################################
# XinLiang Hu
# 2023-08-23

#设置工作目录
setwd("C:/Users/huxin/Desktop/Final_result_analysys/01_Bacteria_species_abundance/Alphan-diversity/")

#加载所需要的R包
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)

#读取数据表
abund_table = read.table("total-7_abundance_table_species_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
metadata = read.table("total-7_metadata1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
abund_table[1:5, 1:5]
head(metadata)
abund_table <- t(abund_table)

###分别计算Shannon指数和Richness
#Shannon指数
alpha_shannon <- diversity(abund_table, "shannon")
head(alpha_shannon)
write.table(Richness,file="all-7_shannon.txt",sep = "\t", quote = FALSE, row.names = T, col.names = FALSE)

#Richness
Richness <- rowSums(abund_table > 0)
head(Richness)
write.table(Richness,file="all-7_richness.txt",sep = "\t", quote = FALSE, row.names = T, col.names = FALSE)

#myread()函数用于读取和处理数据
myread <- function(filename, colname){
  table_tmp = read.table(file = filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(table_tmp) <- c("Individuals", colname)  #添加列名
  box <- merge(table_tmp, metadata, by = "Individuals") #合并物种shannon指数表和对应的metadata
  return(box)
}

#myDiversity()函数用于绘制箱线图
myDiversity <- function(mydata, xdata, ydata, Type, ylabel){
  if(nrow(metadata) == nrow(mydata)){
    ggplot(data = mydata, aes(x=xdata, y=ydata, fill=Type)) +   #绘制箱线图
      geom_boxplot() +
      theme_classic() + 
      scale_fill_manual(values = c('#CD5B45','#00688B'))+
      theme(legend.position = "none", 
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20), 
            axis.title.y = element_text(size = 20)) + 
      scale_x_discrete(breaks = c("Control", "Obesity"), 
                       labels = c("Control", "Obesity")) + 
      geom_signif(comparisons = list(c("Control", "Obesity")), 
                  test = wilcox.test, 
                  map_signif_level = T, 
                  color = "black", textsize = 6) + 
      labs(x = "", y = ylabel)
  } else {
    print("----------- ERROR -------------")
  }
}

shannon_table <- myread(filename = "all-7_shannon.txt", colname = "Shannon")
shannon_draw <- myDiversity(shannon_table, shannon_box$Type, shannon_box$Shannon, shannon_box$Type, ylabel="Shannon Index")
shannon_draw
ggsave(shannon_draw, filename = "all-7_bacteria_Shannon.pdf", width = 5, height = 7)

richness_table <- myread(filename = "all-7_richness.txt", colname = "Richness")
richness_draw <- myDiversity(richness_table, richness_table$Type, richness_table$Richness, shannon_box$Type, ylabel = "Richness Index")
richness_draw
ggsave(richness_draw, filename = "all-7_bacteria_Richness.pdf", width = 5, height = 7)


#不同国家之间的物种多样性比较
ggboxplot(mybox, x = "Country", y = "Shannon", fill = "Type") + 
  scale_fill_manual(values = c('#CD5B45','#00688B'))+
  scale_x_discrete(breaks = c("Denmark", "China", "Sweden", "Spain", "Australia"), 
                   labels = c("Denmark\n(n=278)", "China\n(n=383)", "Sweden\n(n=92)",
                              "Spain\n(n=59)", "Australia\n(n=50)"))+
  theme(legend.position = "top", 
        axis.text.x = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 13)) +
  stat_compare_means(method = "kruskal.test") + 
  labs(x = "", y = "Shannon Index")

ggsave(filename = "all-7_shannon_Country_0.png", width = 7, height = 7, units = "in", dpi = 300)
ggsave(filename = "all-7_shannon_Country_0.pdf", width = 7, height = 7)

#统计检验和P值矫正
with(data = mybox, pairwise.wilcox.test(x=Shannon, g=Country, p.adjust.method = "bonferroni"))


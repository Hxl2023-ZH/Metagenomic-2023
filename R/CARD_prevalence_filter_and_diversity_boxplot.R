
################################################################################
###                            抗性基因分析
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/11_Bacterial_CARD/")

#################################### 过滤数据 ################################
library(dplyr)

df <- read.delim("all-7_CARD_Gene_relative_abundance.txt", row.names = 1, check.names = FALSE)

df[1:5, 1:5]

df1 <- df

#转为0和1的二进制表
df1[df1>0] <- 1
df1[5:10, 859:862]
#write.table(df1,file = "all-7_CARD_binary.txt",sep = "\t",quote = FALSE)

#添加求和列
df2 <- mutate(df,rowsum=apply(df1, 1, sum))
df2[1:5, 1:5]
nrow(df2)

#根据求和列过滤数据，rowsum > 2
df3 <- subset(df2, df2$rowsum > 2)
nrow(df3)
#去除rowsum列，保存文件
df4 <- df3[,-863]
colnames(df4)
ncol(df4)

write.table(df4,file = "all-7_CARD_prevalence_filter.txt",sep = "\t",quote = FALSE)

################################ 计算多样性指数 #######################################
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)

#读取数据表
abund_table = read.delim("all-7_CARD_prevalence_filter.txt", row.names = 1, check.names = FALSE)
metadata = read.table("metadata_sex_age_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
head(metadata)
abund_table <- t(abund_table)

#计算Shannon指数
alpha_shannon <- diversity(abund_table, "shannon")
alpha_simpson <- diversity(abund_table, "simpson")
alpha_richness <- rowSums(abund_table > 0)

#保存文件
write.table(alpha_shannon,file="all-7_CARD_gene_shannon.txt",sep = "\t", quote = FALSE, row.names = T, col.names = FALSE)
write.table(alpha_simpson,file="all-7_CARD_gene_simpson.txt",sep = "\t", quote = FALSE, row.names = T, col.names = FALSE)
write.table(alpha_richness,file="all-7_CARD_gene_richness.txt",sep = "\t", quote = FALSE, row.names = T, col.names = FALSE)


shannon_table = read.table("all-7_CARD_gene_shannon.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
simpson_table = read.table("all-7_CARD_gene_simpson.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
richness_table = read.table("all-7_CARD_gene_richness.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)


colnames(shannon_table) <- c("Individual", "Shannon") ##添加列名
colnames(simpson_table) <- c("Individual", "Simpson")
colnames(richness_table) <- c("Individual", "Richness")

nrow(metadata) == nrow(shannon_table) #查看两个表的行数是否一致

shannon_box <- merge(shannon_table, metadata, by = "Individual") #合并物种shannon指数表和对应的metadata
simpson_box <- merge(simpson_table, metadata, by = "Individual")
richness_box <- merge(richness_table, metadata, by = "Individual")
head(shannon_box)

################################# 绘图 ############################

#准备绘图函数myDraw()
myDraw <- function(data, xdata, ydata, fill=study_condition, ylabe){
  
  ggplot(data, aes(x=xdata, y=ydata, fill=study_condition)) + 
    geom_boxplot() +
    theme_classic() + 
    scale_fill_manual(values = c("#38BBB5","#DCC03D")) +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20)) + 
    scale_x_discrete(breaks = c("Control", "Obesity"), 
                     labels = c("Control", "Obesity")) + 
    geom_signif(comparisons = list(c("Control", "Obesity")), 
                test = wilcox.test, 
                map_signif_level = F, 
                color = "black", textsize = 6) +
    labs(x = "", y = ylabe)
}

shannon_boxplot <- myDraw(shannon_box, shannon_box$study_condition, shannon_box$Shannon, ylabe = "CARD Shannon Index")
simpson_boxplot <- myDraw(simpson_box, simpson_box$study_condition, simpson_box$Simpson, ylabe = "CARD Simpson Index")
richness_boxplot <- myDraw(richness_box, richness_box$study_condition, richness_box$Richness, ylabe = "CARD Richness Index")

shannon_mean <- aggregate(Shannon ~ study_condition, data = mybox, median)
shannon_mean
ggsave(shannon_boxplot, filename = "all-7_CARD_Gene_Shannon.png", width = 5, height = 7, units = "in", dpi = 300)
ggsave(shannon_boxplot, filename = "all-7_CARD_Gene_Shannon.pdf", width = 5, height = 7)
ggsave(simpson_boxplot, filename = "all-7_CARD_Gene_Simpson.pdf", width = 5, height = 7)
ggsave(richness_boxplot, filename = "all-7_CARD_Gene_Richness.pdf", width = 5, height = 7)


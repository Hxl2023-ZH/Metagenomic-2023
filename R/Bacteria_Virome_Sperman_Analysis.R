
################################################################################
###                   计算bacteriome（包含古菌等）和virome之间的相关性
################################################################################
# XinLiang Hu
# 2023-08-23

############################# 根据prevalence过滤原始数据 #######################
bacteriome = read.table("total-7_species_lefse_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
bacteriome[1:5, 1:5]
bacteriome1 <- bacteriome

bacteriome1[bacteriome1>0] <- 1 #转为0和1的二进制表

bacteriome2 <- mutate(bacteriome,rowsum=apply(bacteriome1, 1, sum)) #添加求和列
bacteriome2[1:5, 1:5]
nrow(bacteriome2)

bacteriome3 <- subset(bacteriome2, bacteriome2$rowsum > 173) #根据求和列过滤数据，大于20%的流行率（prevalence）
nrow(bacteriome3)

bacteriome4 <- bacteriome3[,-863] #去除rowsum列，样本总数为862
colnames(bacteriome4)
ncol(bacteriome4)

write.table(bacteriome4, file="bacteria_species_abundance_filter_20.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)


################################################################################
############################# 加载工作数据 #####################################
setwd("C:/Users/huxin/Desktop/Virome_Bacteris_Cor/family_species/Control/")

#加载R包
library(dplyr)
library(Hmisc)
library(data.table)


bacteriome = read.table("Obesity_bacteria_species_abundance.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
bacteriome[1:5, 1:5]

virome = read.table("phage_Obesity_vOTU_family_relative_abundance.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
virome[1:5, 1:5]


############# 计算spearman相关系数 ###########

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
cor1 <- rcorr(t(virome), t(bacteriome), type = 'spearman')
head(cor1)
#阈值筛选
r <- cor1$r
#r[abs(r) < 0.5] <- 0

#筛选 p<0.05
p <- cor1$P
p <- p.adjust(p, method = 'BH') #p值校正
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
z1 <- z

z1 <- as.data.frame(z)
z1[is.na(z1)] <- 0

write.table(data.frame(z1, check.names = FALSE), 'Obesity_phage_family_species_Cor.txt', col.names = NA, sep = '\t', quote = FALSE)

draw1 = read.table("Control_phage_family_species_Cor.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)

head(draw1)

library(pheatmap)

pdf("Control_test.pdf",width = 7,height = 12)
pheatmap(draw1, 
         border_color = "white",
         display_numbers = matrix(ifelse(abs(draw1) > 0.2, "+", ""),nrow = nrow(draw1)),
         number_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_number = 6,
         main = "Bacteria Species And Phage Family Cor In Control",
         #show_rownames = FALSE,
         fontsize_row = 5,
         fontsize_col = 10
         #angle_col = 45
         )

dev.off()

############  进行长宽变换
library(tidyr)

df = read.csv("phage_family_species_Cor.txt", header = TRUE, sep = "\t", check.names = F)
head(df)

df1 <- gather(df, key = "Virome", value = "Cor", -'Species')
head(df1)

write.table(df1, file="Control_phage_family_species_Cor_reshape.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)


#################################### 统计检验 ##################################
#total Cor number
data <- matrix(c(993, 1838, 1130, 1701), nrow = 2)
colnames(data) <- c("Obesity", "Control")
rownames(data) <- c("Correlative", "non-Correlative")
data

fisher.test(data)

#positive and negtive cor number
data1 <- matrix(c(557, 436, 652, 478), nrow = 2)
colnames(data1) <- c("Obesity", "Control")
rownames(data1) <- c("positive-Correlative", "negative-Correlative")
data1

fisher.test(data1)

#对不同的phage family 进行分组相关性统计检验
setwd("C:/Users/huxin/Desktop/Virome_Bacteris_Cor/family_species/")

library(dplyr)
library(tidyr)
library(rstatix)

phage_family <- sample1 <- read.delim('Phage_family_Species_merge.txt', stringsAsFactors = FALSE, check.names = FALSE)
head(phage_family)

phage_statistic <- phage_family %>% group_by(Virome) %>% wilcox_test(Cor ~ Group, detailed = T) # 分组统计检验
phage_statistic_significant <- phage_statistic %>% filter(p<0.05) # 根据P值过滤
phage_statistic_significant

write.table(phage_statistic, file="phage_family_statistics.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)
write.table(phage_statistic_significant, file="phage_family_statistics_significant.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)


species_statistic <- phage_family %>% group_by(Species) %>% wilcox_test(Cor ~ Group, detailed = T) # 分组统计检验
species_statistic_significant <- species_statistic %>% filter(p<0.05) # 根据P值过滤
species_statistic_significant

write.table(species_statistic, file="Species_statistics.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)
write.table(species_statistic_significant, file="Spesies_statistics_significant.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)




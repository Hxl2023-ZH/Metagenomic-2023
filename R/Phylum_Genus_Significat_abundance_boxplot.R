
################################################################################
###                       Phylum和Genus丰度分析
################################################################################
# XinLiang Hu
# 2023-08-23


setwd("C:/Users/huxin/Desktop/")

####################### 提取出具有显著相对丰度差异的Genus #######################
library(tidyr)

df = read.csv("merged_abundance_table_genus.txt", header = TRUE, sep = "\t", check.names = F)
df[1:5, 1:5]
metadata = read.csv("total-7_metadata1.txt", header = TRUE, sep = "\t", check.names = F)
head(metadata)

df1 <- gather(df, key = "Individuals", value = "Abundance", -'Species') #长宽变换
head(df1)

df2 <- merge(df1, metadata, by="Individuals", sort = TRUE) #添加分组情况
head(df2)
nrow(df1)==nrow(df2)

###### 批量统计检验
library(rstatix)

df_statistic <- df2 %>% group_by(Species) %>% wilcox_test(Abundance~ Group, detailed = T) # 分组统计检验
df_statistic_significant <- df_statistic %>% filter(p<0.001) # 根据P值过滤
df_statistic_significant

write.table(df_statistic, file="Genus_relative_abundance_statistics.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)
write.table(df_statistic_significant, file="Genus_relative_abundance_statistics_significant.txt",sep = "\t",row.names=TRUE, col.names=T, quote = FALSE)

###切换至rstudio Terminal(linux终端)，通过以下命令提取提取Genus相对丰度
#cat Genus_relative_abundance_statistics_significant.txt | cut -f2 > Significat_Genus_names.txt
#sed -i 's/estimate/Species/' Significat_Genus_names.txt
#cat Significat_Genus_names.txt | while read LINE; do  sed -n '/'$LINE'/ p' merged_abundance_table_genus.txt >> Significat_Genus_abundance.txt ; done

#### 绘图
Df = read.csv("Significat_Genus_abundance.txt", header = TRUE, sep = "\t", check.names = F)
Df[1:5, 1:5]

Df1 <- gather(Df, key = "Individuals", value = "Abundance", -'Species') #长宽变换
head(Df1)

Df2 <- merge(Df1, metadata, by="Individuals", sort = TRUE) #添加分组情况
head(Df2)
nrow(Df1)==nrow(Df2)

library(ggplot2)
library(ggpubr)

ggplot(Df2, aes(x=Species, y=log10(Abundance), fill=Group)) + 
  geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values = c("#CD5B45","#00688B")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size=15)) + 
  labs(x = "Significant Different Genus", y = "Rlative abundance(log10)")+
  rotate_x_text(90)

ggsave(filename = "Significant_difference_Genus_boxplot.pdf", width = 9, height = 7)













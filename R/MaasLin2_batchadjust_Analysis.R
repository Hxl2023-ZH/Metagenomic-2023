
################################################################################
###                          差异富集物种分析
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/Bacteria_difference_species")

library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
library(Maaslin2)

###load data
mydata = read.table("total-7_abundance_table_species_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
mymeta = read.table("metadata_sex_age_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE, check.names = F)

mydata <- t(mydata)

#直接处理原始数据遇到问题"Feature table does not appear to be either proportions or counts!"，需要进行如下转换
mydata = mydata %>% apply(2, function(x) x / 100)

mydata[1:5, 1:5]
head(mymeta)
table(mymeta$studyID)

#string变量转化为factor
#mymeta[, 'studyID'] <- as.factor(mymeta[, 'studyID'])

###Performing batch (study) effect adjustment
fit_adjust_batch <- adjust_batch(feature_abd = mydata,
                                 batch = "studyID",
                                 covariates = "study_condition",
                                 data = mymeta,
                                 control = list(verbose = FALSE))

mydata_adj <- fit_adjust_batch$feature_abd_adj

#evaluate the effect of batch adjustment

D_before <- vegdist(t(mydata), method = "bray")
D_after <- vegdist(t(mydata_adj), method = "bray")

set.seed(2023)
fit_adonis_before0 <- adonis2(D_before ~ studyID, data = mymeta)


fit_adonis_after <- adonis2(D_after ~ studyID, data = mymeta)


print(fit_adonis_before0$aov.tab)
print(fit_adonis_before1$aov.tab)
print(fit_adonis_before2$aov.tab)
print(fit_adonis_before3$aov.tab)

print(fit_adonis_after$aov.tab)


##################################### 运行maaslin2 #########################################

df_input_data <- t(mydata_adj) #样本为行
fit_data_random = Maaslin2(
  input_data = df_input_data, 
  input_metadata = mymeta, 
  min_prevalence = 0.1,
  output = "all_output_random", 
  fixed_effects = c("study_condition", "MF", "AgeGroup"),plot_heatmap = TRUE,
  reference = c(("study_condition,Obesity"),("Country,China"),("AgeGroup,Age1"))
  )

#提取和过滤结果
fit_data_df <- as.data.frame(fit_data_random$results)

fit_data_df1 <- fit_data_df %>% filter(metadata == "study_condition") %>% 
  filter(pval < 0.01) %>% 
  filter(qval < 0.01) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature))

write.table(fit_data_df1,file="sin_difference_species_0.01.txt",sep = "\t", quote = FALSE, row.names = T, col.names = TRUE)

df <- read.delim("sin_difference_species_0.01.txt")
df <- df %>% mutate(feature = factor(feature, levels = feature))
df$Group <- as.factor(df$Group)
head(df)
#使用ggplot绘制条形图
ggplot(df, aes(coef, feature, fill=Group)) + 
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("#2C5094","#F45268"))+
  theme_bw() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        legend.position = c(0.92, 0.05),
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25))+
  labs(x = "Coef", y = "")

ggsave(filename = "sin_difference_species.png", width = 20, height = 20, units = "in", dpi = 300)

ggsave(filename = "sin_difference_species1.pdf", width = 20, height = 20)



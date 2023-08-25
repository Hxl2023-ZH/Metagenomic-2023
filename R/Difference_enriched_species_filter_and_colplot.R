
################################################################################
###                   根据prevalence和mean relative abundance过滤差异物种
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/")


df <- read.delim('sin_species_abundance.txt', row.name = 1, check.names = FALSE)
df[1:5,1:5]
rownames(df)

#计算每行feature的平均相对丰度
df_average <- apply(df, 1, FUN=mean)
head(df_average)

write.table(df_average,file="sin_species_average_abundance.txt",sep = "\t", quote = FALSE, row.names = T, col.names = TRUE)

#生成二进制table
df1 <- df
df1[df1 > 0] <- 1
df1[1:5,1:5]

#计算prevarence
df1_sum <- apply(df1, 1, FUN=sum)
df1_sum

write.table(df1_sum,file="sin_species_prevalence.txt",sep = "\t", quote = FALSE, row.names = T, col.names = TRUE)

#在桌面上修改文件夹后重新导入
df2 <- read.delim('sin_species_average_abundance.txt', check.names = FALSE)
head(df2)
df3 <- read.delim('sin_species_prevalence.txt', check.names = FALSE)
head(df3)
df4 <- read.delim('sin_difference_species_0.01.txt', check.names = FALSE)
head(df4)

#合并数据
Df1 <- merge(df2, df3, by="feature")
head(Df1)
Df2 <- merge(Df1, df4, by="feature")
head(Df2)
write.table(Df2,file="sin_difference_species_001_plus.txt",sep = "\t", quote = FALSE, row.names = T, col.names = TRUE)

#过滤数据，prevarence>20%
library(tidyverse)

Df3 <- Df2 %>% filter(metadata == "study_condition") %>% 
  filter(prevarence > 173) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature))

write.table(Df3,file="sin_difference_species_001_filter.txt",sep = "\t", quote = FALSE, row.names = T, col.names = TRUE)

#df <- read.delim("sin_difference_species_0.01.txt")
Df3 <- Df3 %>% mutate(feature = factor(feature, levels = feature))
Df3$Group <- as.factor(Df3$Group)
head(Df3)
#使用ggplot绘制条形图
ggplot(Df3, aes(coef, feature, fill=Group)) + 
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

ggsave(filename = "sin_difference_species_filter.png", width = 20, height = 20, units = "in", dpi = 300)

ggsave(filename = "sin_difference_species_filter.pdf", width = 20, height = 20)

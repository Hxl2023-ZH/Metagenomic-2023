
################################################################################
###                        绘制差异富集物种条形图
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/Bacteria_difference_species_02/Lefse_Maasslin2/")

library(ggplot2)
library(ggpubr)

df <- read.delim("lefse_maasslin2_plot.txt", check.names = F)
head(df)
#df <- df %>% mutate(Group = factor(Group, levels = c("Obesity", "Control")))
#df$Group <- as.factor(df$Group, levels=c("Obesity", "Control"))
df <- df %>% mutate(Features = factor(Features, levels = Features))
df$Group <- as.factor(df$Group)
head(df)
####使用ggplot绘制条形图
ggplot(df, aes(x=COEF, y=Features, fill=Group)) + #构建ggplot绘图几何对象
  geom_col(width = 0.8) +#绘制条形图并设置大小
  scale_fill_manual(values = c('#CD5B45','#00688B'))+#手动设置颜色
  theme_minimal()+#设置主题
  theme(panel.grid =element_blank())+#去除网格线
  #theme(panel.border = element_blank())+#去除边框
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        #legend.position = c( 0.03, 0.97),
        legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16))+
  geom_text(data=subset(df, COEF<0), #设置正值的注释文本的位置和大小
            aes(x=0, y= Features,label= paste0(' ', Features)), 
            size = 6, hjust = 'outward' ) +
  geom_text(data=subset(df, COEF>0), #设置负值的注释文本的位置和大小
            aes(x=0, y= Features,label= paste0(Features, ' ')), 
            size = 6, hjust = 'inward' ) +
  geom_vline(xintercept = c(-2, -1, 0, 1, 2, 2),lty="dashed")+#添加垂直线
  labs(x = "COEF", y = "")#设置X,Y轴标题文本

ggsave(filename = "Maasslin2+lefse_species.png", width = 12, height = 16, units = "in", dpi = 300)
ggsave(filename = "Maasslin2+lefse_species.pdf", width = 12, height = 16)


#hxl
#含抗性基因vOTU相对丰度箱线图的绘制
setwd("C:/Users/huxin/Desktop/Virome_resfinder/")

library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(tidyverse)

dat <- read.delim('res_vOTU_abundance.txt')# row: samples
head(dat)

#数据格式转为长格式，tidyverse
pdat <- dat %>%
  pivot_longer(cols = !c(Individuals,Group),
               names_to = "amrs",
               values_to = "count")
head(pdat)

#绘制箱线图
draw1 <- ggplot(pdat, aes(x=reorder(amrs, count),y=count))+
  geom_col(aes(fill = Group), width = 0.5, position = "dodge") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        legend.position = 'top',
        axis.title.x = element_text(size = 5),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 10)) +
  stat_compare_means(aes(group = Group),method = "wilcox.test",
                     label = "p.signif",hide.ns = TRUE)+
  labs(x = "", y = 'Relative Abundance\n')+
  rotate_x_text(45)

draw1
#保存结果

ggsave(draw1, filename = "vOTU_resfinder_col_dodge.png", width = 12, height = 7, units = "in", dpi = 300)


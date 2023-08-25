
################################################################################
###                   拟杆菌门和厚壁菌门分析
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/01_Bacteria_species_abundance/High_Classic")

library(wesanderson)
library(colortools)
library(ggpubr)
library(ggsignif)
library(ggplot2)
library(patchwork)

# 导入数据并查看数据集格式
Class <- read.delim('Bacteroidetes_Firmcutes2.txt', row.name = 1, check.names = FALSE)
metadata <- read.delim('total-7_metadata1.txt', check.names = FALSE)
head(metadata)
head(Class)
#绘图
p1 <- ggplot(Class, aes(Group,log10(F_B))) +
  geom_boxplot(aes(fill = Group), width = 0.6) +
  scale_fill_manual(values = c('#CD5B45','#00688B'))+
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        legend.position = 'none',
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  geom_signif(comparisons = list(c("Control", "Obesity")), 
              test = wilcox.test, 
              map_signif_level = T, 
              color = "black", textsize = 6) +
  labs(x = "", y = 'Firmicutes/Bacteroidetes Ratio(log10)')
p1

p2 <- ggplot(Class, aes(Group,Firmicutes)) +
  geom_boxplot(aes(fill = Group), width = 0.6) +
  scale_fill_manual(values = c('#CD5B45','#00688B'))+
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        legend.position = 'none',
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  geom_signif(comparisons = list(c("Control", "Obesity")), 
              test = wilcox.test, 
              map_signif_level = T, 
              color = "black", textsize = 6) +
  labs(x = "", y = 'Firmicutes Relative Abundance')
p2
p3 <- ggplot(Class, aes(Group,Bacteroidetes)) +
  geom_boxplot(aes(fill = Group), width = 0.6) +
  scale_fill_manual(values = c('#CD5B45','#00688B'))+
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        legend.position = 'none',
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15)) +
  geom_signif(comparisons = list(c("Control", "Obesity")), 
              test = wilcox.test, 
              map_signif_level = T, 
              color = "black", textsize = 6) +
  labs(x = "", y = 'Bacteroidetes Relative Abundance')
p3

p <- p1+p2+p3
p
ggsave(p, filename = "Bacteroidetes-Firmicutes.png", width = 7, height = 7, units = "in", dpi = 300)

ggsave(p, filename = "Bacteroidetes-Firmicutes.pdf", width = 7, height = 7)



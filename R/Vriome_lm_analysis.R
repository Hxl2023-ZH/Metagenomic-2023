
################################################################################
###                       vOTU和物种多样性的相关分析
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/09_Virome/Virome_workspace/Virome_lm_analysis")

library(ggpmisc)
library(ggplot2)
library(ggpubr)

df = read.table("virome_bacteria_shannon_lm.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
head(df)

#Remove confidence region (se = FALSE)
# Extend the regression lines: fullrange = TRUE
draw1 <- ggplot(df, aes(x = Bactriome_Shannon, y = vOTU_Shannon))+
  geom_point(aes(color = Group, shape = Group)) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))+
  #geom_rug(aes(color =Group)) +
  geom_smooth(aes(color = Group),
              method = lm, 
              se = TRUE, 
              fullrange = TRUE)+
  scale_color_manual(values = c('#CD5B45','#00688B'))+
  ggpubr::stat_cor(method = "spearman", aes(color = Group), label.x = 1)+
  labs(x = "Shannon index (Bacteriome)\n", y = 'Shannon index (Virome)\n')
  
draw1
ggsave(draw1, filename = "virome_bacteriome_lm.png", width = 7, height = 7, units = "in", dpi = 300)
ggsave(draw1, filename = "virome_bacteriome_lm.pdf", width = 7, height = 7)



 

################################################################################
###                                噬菌体基因组特征分析
################################################################################
# XinLiang Hu
# 2023-08-23


setwd("C:/Users/huxin/Desktop/Final_result_analysys/09_Virome/Virome_workspace/vOTU_Phage/Family/Phage_GC_Size/")

library(ggplot2)

############### 绘图 #############

library(ggplot2)
library(ggsignif)
library(ggpubr)

head(df)
#length
ggplot(df, aes(x=phylum, y=Length_kbp, fill=phylum)) + 
  geom_boxplot() +
  geom_jitter(size = 0.2, alpha = 0.05)+
  theme_classic() + 
  #scale_fill_manual(values = c("#2C797B", "#428273", "#588B6B", "#6E9364", "#849C5C", "#9AA554")) +
  scale_fill_manual(values = c('#00688B', '#588B6B', '#B28F7F', '#FFA07A', '#CD5B45', '#FFD700')) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20)) + 
  geom_signif(comparisons = list(c("Actinobacteria", "Bacillota"),
                                 c("Actinobacteria", "Bacteroidota"),
                                 c("Actinobacteria","Firmicutes"),
                                 c("Actinobacteria","Proteobacteria"),
                                 c("Actinobacteria", "Other"),
                                 c("Bacteroidota", "Bacillota"),
                                 c("Bacteroidota","Firmicutes"),
                                 c("Bacteroidota","Proteobacteria"),
                                 c("Bacteroidota", "Other")),
              y_position = c(260,250,240,220,230,210,202,186,194),
              test = wilcox.test,
              map_signif_level = T, 
              color = "black",
              textsize = 3.5)+
  labs(x = "Host", y = "Phage size Content(kbp)")+
  rotate_x_text(90)

ggsave(filename = "Host_Phage_size1.png", width = 9, height = 7, units = "in", dpi = 300)
ggsave(filename = "Host_Phage_size1.pdf", width = 9, height = 7)

#GC
ggplot(df, aes(x=phylum, y=GC, fill=phylum)) + 
  geom_boxplot() +
  geom_jitter(size = 0.2, alpha = 0.05)+
  theme_classic() + 
  #scale_fill_manual(values = c("#2C797B", "#428273", "#588B6B", "#6E9364", "#849C5C", "#9AA554")) +
  scale_fill_manual(values = c('#00688B', '#588B6B', '#B28F7F', '#FFA07A', '#CD5B45', '#FFD700')) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20)) + 
  geom_signif(comparisons = list(c("Actinobacteria", "Bacillota"),
                                 c("Actinobacteria", "Bacteroidota"),
                                 c("Actinobacteria","Firmicutes"),
                                 c("Actinobacteria","Proteobacteria"),
                                 c("Actinobacteria", "Other"),
                                 c("Bacteroidota", "Bacillota"),
                                 c("Bacteroidota","Firmicutes"),
                                 c("Bacteroidota","Proteobacteria"),
                                 c("Bacteroidota", "Other")),
              y_position = c(80,78,76,74,72,70,68,64,66),
              test = wilcox.test,
              map_signif_level = T, 
              color = "black",
              textsize = 3.5)+
  labs(x = "Host", y = "Phage GC Content")+
  rotate_x_text(90)

ggsave(filename = "Host_Phage_GC1.png", width = 9, height = 7, units = "in", dpi = 300)
ggsave(filename = "Host_Phage_GC1.pdf", width = 9, height = 7)





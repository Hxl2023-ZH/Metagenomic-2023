
################################################################################
###                        不同物种分类水平柱状图绘制
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/01_Bacteria_species_abundance/High_Classic/")

# 导入数据并查看数据集格式
Class <- read.delim('merged_abundance_table_phy.txt', row.name = 1, check.names = FALSE)
metadata <- read.delim('total-7_metadata1.txt', check.names = FALSE)
head(metadata)

Class[1:5, 1:5]
# 相对丰度
#Order <- as.data.frame(lapply(Order, function(x) x / 100))
#row.names(phylum_per) <- row.names(phylum) #加一下行名
# 计算每个门水平的平均丰度 便于后续筛选                                 
Class.ave <- apply(Class, 1, FUN=mean)
Class.2 <- cbind(Class,Class.ave)[order(-Class.ave),] #排个序
Class.2[, 1]
# 选择丰度最高的10个门 剩下的放入others里
Class.2 <- subset(Class.2[1:10,], select=-Class.ave)
# 统计others丰度
Class.2 <- rbind(Class.2, Others=apply(Class.2, 2, function(x){100-sum(x)}))
# 加一列行名 便于后续的长宽转换
Class.2 <- cbind(ClassID=row.names(Class.2), Class.2)
Class.2[1:5, 1:5]

# 长宽转换
library(reshape2)
# 因子排个序
Class.2$ClassID <- factor(Class.2$ClassID, levels = rev(Class.2$ClassID))
Class.gg <- melt(Class.2, id.vars="ClassID", variable.name="SampleID", value.name="Abundance")
head(Class.gg)

#添加分组
df <- merge(Class.gg, metadata, by="SampleID")
head(df)


library(wesanderson)
library(colortools)
library(ggpubr)
library(ggsignif)

ggplot(df, aes(x=SampleID,y=Abundance))+
  geom_col(aes(fill = ClassID), width = 0.8, position = "stack") +
  facet_grid(~ Group, scales = "free_x", space='free') + 
  scale_fill_manual(values=c("gray",rev(wheel("red")[1:10]))) +
  #scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        legend.position = 'right',
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10)) +
  guides(fill=guide_legend(title="Top Phylum"))+
  labs(x = "Individuals", y = 'Relative Abundance')

ggsave(filename = "Genus_relative_abundance.png", width = , height = 9, units = "in", dpi = 300)
ggsave(filename = "Bacteroides_vulgatus_57955_SNPs_frequence.pdf", width = 5, height = 7)








################################################################################
###                                PCoA分析
################################################################################
# XinLiang Hu
# 2023-08-23


setwd("C:/Users/huxin/Desktop/pcoa_work/species_pcoa/")

#加载数据，丰度表格行为样本名，列为物种名
otu <- read.delim('total-7_abundance_table_species_0.txt', row.names = 1,check.names=F)
otu[1:3, 1:3]
metadata <- read.delim('metadata_sex_age_0.txt', stringsAsFactors = FALSE,check.names=F)
head(metadata)
library(vegan)
library(ggplot2)

#计算Brar_curtis
otu.distance <- vegdist(otu)

#pcoa分析
pcoa <- cmdscale (otu.distance,eig=TRUE) #对距离进行排序
pc12 <- pcoa$points[,1:2] #提取样本点前两坐标轴
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2) #计算前两坐标轴解释量

#数据格式转换及数据整合
pc12 <- as.data.frame(pc12)
pc12$Individuals <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p
#colnames(group) <- c("samples","group","Country")
colnames(metadata)

df <- merge(pc12,metadata,by="Individuals")
head(df)
color=c("#1597A5","#FFC24B",'#CD5B45',"#100D01",'#34E10E')
p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=Country,shape=Country))+
  theme_bw()+
  geom_point(size=1.8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PCoA1 ",pc[1],"%"),
       y=paste0("PCoA2 ",pc[2],"%"))+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#1597A5","#FFC24B",'#CD5B45',"#100D01",'#34E10E'))+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.position = c(0.9, 0.85),
        panel.grid=element_blank())+ 
  stat_ellipse(data=df,geom = "polygon",level=0.95,linetype = 2,size=0.5,aes(fill=Country),alpha=0.2,show.legend = T)
p1
#按分组进行PERMANOVA分析
out.adoins0 = adonis2(otu.distance~df$V2, distance = 'bray', permutations = 999)
out.adoins1 = adonis2(otu ~ Group, data = metadata, permutations = 999, distance = 'bray')
out.adoins2 = adonis2(otu ~ Country, data = metadata, permutations = 999, distance = 'bray')
out.adoins3 = adonis2(otu ~  studyID, data = metadata, permutations = 999, distance = 'bray')
out.adoins4 = adonis2(otu ~ MF, data = metadata, permutations = 999, distance = 'bray')
out.adoins5 = adonis2(otu ~ AgeGroup, data = metadata, permutations = 999, distance = 'bray')
out.adoins0
out.adoins1
out.adoins2
out.adoins3
out.adoins4
out.adoins5

p1 + annotate('text', label='PERMANOVA:\n***', x=0.57, y=0.21, size=3.5)
ggsave(filename = "Bacteria_Country_PCoA.png", width = 7, height = 7, units = "in", dpi = 300)




################################################################################
###                                细菌毒力因子逻辑回归分析
################################################################################
# XinLiang Hu
# 2023-08-23

setwd("C:/Users/huxin/Desktop/Final_result_analysys/10_ARG_Vfs_glm")
library(dplyr)

df <- read.table("tmp.txt", header=T, check.names = FALSE)
head(df)


#df$Group<-factor(df$Group, levels = c(1,0), labels = c("Obesity","Control")) #转换数据类型
df$Group1<-factor(df$Group1, levels = c(0,1), labels = c("Control","Obesity")) #转换数据类型

#df$Group1<-factor(df$Group1)
#df$Country<-factor(df$Country)
df$studyID <- factor(df$studyID)
df$AgeGroup <- factor(df$AgeGroup)
str(df)

fit.full<-glm(Group1~Bacterial_Vfs_Shannon+Bacterial_Vfs_Simpson+Bacterial_Vfs_Richness+MF+
                studyID+AgeGroup, 
              data=df, family = binomial()) #拟合逻辑回归分析

fit.result<-summary(fit.full)
fit.result

fit.full1 <- step(fit.full, direction = "both") #逐步回归分析
#fit.full1 <- step(fit.full, direction = "both", scope = list(lower=Group~Bacterial_Vfs_Shannon+Bacterial_Vfs_Richness, upper=fit.full)) #逐步回归分析

fit.result1<-summary(fit.full1)
fit.result1


df1<-exp(coefficients(fit.result1)) #对结果进行指数化
df2 <- exp(confint(fit.full1))
df3<-cbind(df1,df2)
df4<-data.frame(df3[-1,c(1,4,5,6)])
df4$Var<-rownames(df4)
colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")
df5<-df4[,c(5,1,2,3,4)]
df5$OR_mean<-df5$OR
df5$OR<-paste0(round(df5$OR,4),
               "(",
               round(df5$OR_1,4),
               "~",
               round(df5$OR_2,4),
               ")")
df5$Pvalue <- log(df5$Pvalue)
df5$Pvalue<-round(df5$Pvalue,3)
df5

write.csv(df5,file = "Arg_VFs_forestplot_data_0904.csv",
          quote = F,row.names = F)

#绘制森林图，可跳过
library(forestplot)

#需先调整表格
fp<-read.csv("Arg_VFs_forestplot_data.csv",header=T)

forestplot(labeltext=as.matrix(fp[,1:3]),
           mean=fp$OR_mean,
           lower=fp$OR_1,
           upper=fp$OR_2,
           zero=0,
           boxsize=0.2,
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,
           lwd.ci=2, 
           col=fpColors(box='#C25E5E',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="Odds Ratio",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = "solid",
           title = "Forest Plot", 
           line.margin = 0.08,
           graph.pos=2)




################################################################################
###                                细菌毒力因子回归分析
################################################################################
# XinLiang Hu
# 2023-08-23
#以下内容主要参考https://zhuanlan.zhihu.com/p/344796817

setwd("C:/Users/huxin/Desktop/Final_result_analysys/ARG_Vfs_glm/")
library(dplyr)

df <- read.table("all-7_Arg_VFs_for_glm.txt", header=T, check.names = FALSE)
head(df)


df$Group<-factor(df$Group, levels = c(0,1), labels = c("Obesity","Control"))
head(df)

fit.full<-glm(Group~Bacterial_Vfs_Shannon+Bacterial_Vfs_Simpson+
                Country+MF+AgeGroup, 
              data=df, family = binomial())

fit.result<-summary(fit.full)
fit.result

df1<-fit.result$coefficients
df2<-confint(fit.full)
df3<-cbind(df1,df2)
df4<-data.frame(df3[-1,c(1,4,5,6)])
df4$Var<-rownames(df4)
colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")
df5<-df4[,c(5,1,2,3,4)]
df5$OR_mean<-df5$OR
df5$OR<-paste0(round(df5$OR,2),
               "(",
               round(df5$OR_1,2),
               "~",
               round(df5$OR_2,2),
               ")")
df5$Pvalue<-round(df5$Pvalue,3)
df5

write.csv(df5,file = "Arg_VFs_forestplot_data.csv",
          quote = F,row.names = F)


library(forestplot)
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





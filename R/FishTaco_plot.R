
################################################################################
###                        FishTaco分析结果绘图
################################################################################
# XinLiang Hu
# 2023-08-23
#分析所用的代码主要来源于FishTaco官网http://borenstein-lab.github.io/fishtaco/visualization.html

setwd("C:/Users/huxin/Desktop/Final_result_analysys/De_novo_inference/Pathway_Multi_taxa/Pathway_Control_Enriched/")

library(ggplot2)
library(FishTacoPlot)

p <- MultiFunctionTaxaContributionPlots(input_dir="C:/Users/huxin/Desktop/Final_result_analysys/De_novo_inference/Pathway_Multi_taxa/Pathway_Control_Enriched/", 
                                        input_prefix="Pathway_Control",
                                        input_taxa_taxonomy="taxontable.tab", 
                                        sort_by="list", 
                                        plot_type="bars",
                                        #input_function_filter_list=c("ko00020", "ko00540","ko02040"), 
                                        add_predicted_da_markers=TRUE, 
                                        add_original_da_markers=TRUE)
p

p + 
  #scale_x_continuous(breaks=c(1, 2, 3), labels=c("TCA cycle", "LPS\nbiosynthesis", "Flagellar\nassembly")) +
  guides(fill=guide_legend(nrow=7)) + 
  ylab("Wilcoxon test statistic (W)") +
  theme(plot.title=element_blank(), 
        axis.title.x=element_text(size=12,colour="black",face="plain"),
        axis.text.x=element_text(size=10,colour="black",face="plain"), 
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=10,colour="black",face="plain"), 
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major.x = element_line(colour="light gray"), 
        panel.grid.major.y = element_line(colour="light gray"),
        panel.grid.minor.x = element_line(colour="light gray"), 
        panel.grid.minor.y = element_line(colour="light gray"),
        panel.background = element_rect(fill="transparent",colour=NA), 
        panel.border = element_rect(fill="transparent",colour="black"),
        legend.background=element_rect(colour="black"), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=8,face="plain"),
        legend.key.size=unit(0.8,"line"), 
        #legend.margin=unit(0.1,"line"), 
        legend.position="bottom")


ggsave(filename = "Fishtaco_Control_result.png", width = 14, height = 7, units = "in", dpi = 300)
ggsave(filename = "Fishtaco_Control_result.pdf", width = 14, height = 7)



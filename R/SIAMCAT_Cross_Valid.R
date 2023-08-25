
################################################################################
###                                机器学习
################################################################################
# XinLiang Hu
# 2023-08-23
#代码主要参考SIAMCAT官方文档


setwd("/home/hxl/work_space7/Combined_ROC/Single_Cross_valid")

library(SIAMCAT)
library(dplyr)
library(ggplot2)
#features (in rows) x samples (in columns)
#Please note that SIAMCAT is supposed to work with relative abundances. 
#Other types of data (e.g. counts) will also work, but not all functions of the package will result in meaningful outputs.

mydata = read.table("all-7_kos_mapped_relative_abundance.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
mymeta = read.table("metadata_sex_age_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
head(mymeta)

rownames(mymeta)
rownames(mymeta) <- mymeta$Individual
rownames(mymeta)
#mydata <- t(mydata)
#mydata = mydata %>% apply(2, function(x) x / 100)
#mydata <- as.matrix(mydata)
# check that metadata and features agree
stopifnot(all(colnames(mydata) == mymeta$Individual))
# datasets
datasets <- c('China-2017', 'Australia-2015', 'Denmark-2013', 'China-2020', 'China-2021', 'Sweden-2013')

table(mymeta$studyID, mymeta$study_condition)

###################################### Compare Associations ###############################
#Compute Associations with SIAMCAT
assoc.list <- list()
for (d in datasets){
  # filter metadata and convert to dataframe
  meta.train <- mymeta %>% 
    filter(studyID==d) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Individual
  
  # create SIAMCAT object
  mylabel <- create.label(meta=meta.train, label='study_condition', case='Obesity')
  sc.obj <- siamcat(feat=mydata, label=mylabel, meta=meta.train)
  #sc.obj <- siamcat(feat=mydata, meta=meta.train, label='study_condition', case='Obesity')
  
  # test for associations
  sc.obj <- check.associations(sc.obj, log.n0=1e-06, feature.type = 'original')
  
  # extract the associations and save them in the assoc.list
  temp <- associations(sc.obj)
  temp$genus <- rownames(temp)
  assoc.list[[d]] <- temp %>% 
    select(genus, fc, auc, p.adj) %>% 
    mutate(studyID=d)
}
# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(genus!='unclassified')
head(df.assoc)
#Plot Heatmap for Interesting Genera
genera.of.interest <- df.assoc %>% 
  group_by(genus) %>% 
  summarise(m=mean(auc), n.filt=any(auc < 0.2 | auc > 0.6), 
            .groups='keep') %>% 
  filter(n.filt) %>% 
  arrange(m)

df.assoc %>% 
  # take only genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # convert to factor to enforce an ordering by mean AUC
  mutate(genus=factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  # convert to factor to enforce ordering again
  mutate(studyID=factor(studyID, levels = datasets)) %>% 
  # annotate the cells in the heatmap with stars
  mutate(l=case_when(p.adj < 0.05~'*', TRUE~'')) %>%  
  ggplot(aes(y=genus, x=studyID, fill=fc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
                       limits=c(-2.7, 2.7), name='Generalized\nfold change') + 
  theme_minimal() + 
  geom_text(aes(label=l)) +
  theme(panel.grid = element_blank()) + 
  xlab('') + ylab('') +
  theme(axis.text = element_text(size=7))

############# Study as Confounding Factor
df.meta <- mymeta %>% 
  as.data.frame()
rownames(df.meta) <- df.meta$Individual
sc.obj <- siamcat(feat=mydata, meta=df.meta, label='study_condition', case='Obesity')

check.confounders(sc.obj, fn.plot = './confounder_plot_meta.pdf',
                  feature.type='original')


####################################### ML Meta-analysis ####################################
#########Train LASSO Models
# create tibble to store all the predictions
auroc.all <- tibble(study.train=character(0), 
                    study.test=character(0),
                    AUC=double(0))
# and a list to save the trained SIAMCAT objects
sc.list <- list()
for (i in datasets){
  # restrict to a single study
  meta.train <- mymeta %>% 
    filter(studyID==i) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Individual
  
  ## take into account repeated sampling by including a parameters
  ## in the create.data.split function
  ## For studies with repeated samples, we want to block the
  ## cross validation by the column 'Individual_ID'
  block <- NULL
  
  # create SIAMCAT object
  sc.obj.train <- siamcat(feat=mydata, meta=meta.train, 
                          label='study_condition', case='Obesity')
  # normalize features
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
                                     norm.param=list(log.n0=1e-05, sd.min.q=0),feature.type = 'original')
  # Create data split
  sc.obj.train <- create.data.split(sc.obj.train, num.folds = 10, num.resample = 10, inseparable = block)
  # train LASSO model
  sc.obj.train <- train.model(sc.obj.train, method='lasso')
  
  ## apply trained models to other datasets
  
  # loop through datasets again
  for (i2 in datasets){
    if (i == i2){
      # make and evaluate cross-validation predictions (same dataset)
      sc.obj.train <- make.predictions(sc.obj.train)
      sc.obj.train <- evaluate.predictions(sc.obj.train)
      auroc.all <- auroc.all %>% 
        add_row(study.train=i, study.test=i,
                AUC=eval_data(sc.obj.train)$auroc %>% as.double())
    } else {
      # make and evaluate on the external datasets
      # use meta.ind here, since we want only one sample per subject!
      meta.test <- mymeta %>% 
        filter(studyID==i2) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Individual
      sc.obj.test <- siamcat(feat=mydata, meta=meta.test,
                             label='study_condition', case='Obesity')
      # make holdout predictions
      sc.obj.test <- make.predictions(sc.obj.train, 
                                      siamcat.holdout = sc.obj.test)
      sc.obj.test <- evaluate.predictions(sc.obj.test)
      auroc.all <- auroc.all %>% 
        add_row(study.train=i, study.test=i2,
                AUC=eval_data(sc.obj.test)$auroc %>% as.double())
    }
  }
  # save the trained model
  sc.list[[i]] <- sc.obj.train
}

test.average <- auroc.all %>% 
  filter(study.train!=study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC=mean(AUC), .groups='drop') %>% 
  mutate(study.train="Average")


###Now that we have the AUROC values, we can plot them into a nice heatmap
# combine AUROC values with test average
bind_rows(auroc.all, test.average) %>% 
  # highlight cross validation versus transfer results
  mutate(CV=study.train == study.test) %>%
  # for facetting later
  mutate(split=case_when(study.train=='Average'~'Average', TRUE~'none')) %>% 
  mutate(split=factor(split, levels = c('none', 'Average'))) %>% 
  # convert to factor to enforce ordering
  mutate(study.train=factor(study.train, levels=c(datasets, 'Average'))) %>% 
  mutate(study.test=factor(study.test, 
                           levels=c(rev(datasets),'Average'))) %>% 
  ggplot(aes(y=study.test, x=study.train, fill=AUC, size=CV, color=CV)) +
  geom_tile() + theme_minimal() +
  # text in tiles
  geom_text(aes_string(label="format(AUC, digits=2)"), 
            col='black', size=4)+
  # color scheme
  scale_fill_gradientn(colours=rev(c('darkgreen','forestgreen', 
                                     'chartreuse3','lawngreen', 
                                     'red')), limits=c(0.20, 1)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1, size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('vOTU Training Set\n') + ylab('Test Set\n') + 
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 1), guide=FALSE) + 
  facet_grid(~split, scales = 'free', space = 'free')

ggsave(filename = "KOs_cross_vaild.png", width = 7, height = 7, units = "in", dpi = 300)
ggsave(filename = "KOs_cross_vaild.pdf", width = 7, height = 7)

#################################### Investigate Feature Weights ###########################
#delete SPA-2014
weight.list <- list()
for (d in datasets){
  sc.obj.train <- sc.list[[d]]
  # extract the feature weights out of the SIAMCAT object
  temp <- SIAMCAT::feature_weights(sc.obj.train)
  temp$genus <- rownames(temp)
  # save selected info in the weight.list
  weight.list[[d]] <- temp %>% 
    select(genus, median.rel.weight, mean.rel.weight, percentage) %>% 
    mutate(studyID=d) %>% 
    mutate(r.med=rank(-abs(median.rel.weight)), 
           r.mean=rank(-abs(mean.rel.weight)))
}
# combine all feature weights into a single tibble
df.weights <- bind_rows(weight.list)
df.weights <- df.weights %>% filter(genus!='unclassified')

# compute absolute feature weights
abs.weights <- df.weights %>% 
  group_by(studyID) %>% 
  summarise(sum.median=sum(abs(median.rel.weight)),
            sum.mean=sum(abs(mean.rel.weight)),
            .groups='drop')

df.weights %>% 
  full_join(abs.weights) %>% 
  # normalize by the absolute model size
  mutate(median.rel.weight=median.rel.weight/sum.median) %>% 
  # only include genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # highlight feature rank for the top 20 features
  mutate(r.med=case_when(r.med > 20~NA_real_, TRUE~r.med)) %>%
  # enforce the correct ordering by converting to factors again
  mutate(genus=factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  mutate(studyID=factor(studyID, levels = datasets)) %>% 
  ggplot(aes(y=genus, x=studyID, fill=median.rel.weight)) + 
  geom_tile() + 
  scale_fill_gradientn(colours=rev(
    c('#007A53', '#009F4D', "#6CC24A", 'white',
      "#EFC06E", "#FFA300", '#BE5400')), 
    limits=c(-0.15, 0.15)) +
  theme_minimal() + 
  geom_text(aes(label=r.med), col='black', size= 2) +
  theme(panel.grid = element_blank()) + 
  xlab('') + ylab('') +
  theme(axis.text = element_text(size=6))














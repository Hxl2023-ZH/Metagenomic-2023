
################################################################################
###                                机器学习
################################################################################
# XinLiang Hu
# 2023-08-23
#代码主要参考SIAMCAT官方文档

setwd("C:/Users/huxin/Desktop/Final_result_analysys/Siamcat_lasso_10")
library(SIAMCAT)
library(dplyr)

##################################### 加载数据，构架SIAMCAT所需的对象 #####################
feat.crc.zeller = read.table("all-7_vOTU_relative_abudance.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
meta.crc.zeller = read.delim("metadata_sex_age_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE, check.names = F)
#出现报错：Error in validObject(.Object) : 类别为“phyloseq”的对象不对:
#运行以下命令处理报错
rownames(meta.crc.zeller)
rownames(meta.crc.zeller) <-meta.crc.zeller$Individual
rownames(meta.crc.zeller)

sample_names(feat.crc.zeller)

stopifnot(all(colnames(feat.crc.zeller) == meta.crc.zeller$Individual))

feat.crc.zeller[1:3, 1:3]
dim(feat.crc.zeller)
head(meta.crc.zeller)

#设置标签
label.crc.zeller <- create.label(meta=meta.crc.zeller, label='study_condition', case='Obesity')
#构建对象
sc.obj <- siamcat(feat=feat.crc.zeller, label=label.crc.zeller, meta=meta.crc.zeller)

show(sc.obj)
#过滤数据集
sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 0.00001)
show(sc.obj)

###################################### Association Testing ################################
sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', 
                 panels = c('fc', 'prevalence', 'auroc'))

##################################### Confounder Testing ###################################
check.confounders(sc.obj, fn.plot = 'confounder_plots.pdf',
                  meta.in = NULL, feature.type = 'filtered')


##################################### Model Building #######################################
#Data Normalization
sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

#Prepare Cross-Validation
sc.obj <-  create.data.split(sc.obj, num.folds = 10, num.resample = 10)

#Model Training
sc.obj <- train.model(sc.obj, method = "lasso")

# get information about the model type
model_type(sc.obj)
# access the models
models <- models(sc.obj)
models[[1]]$model

#################################### Make Predictions #######################################
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)

head(pred_matrix)

################################## Model Evaluation and Interpretation #######################
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot('test1'=sc.obj)
model.interpretation.plot(sc.obj, fn.plot = 'interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')





################################################################################################
############################ KEGG abundance table ##############################################

##################################### 加载数据，构架SIAMCAT所需的对象 #####################

feat.crc.zeller1 = read.table("all-7_kos_mapped_relative_abundance.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
#feat.crc.zeller1 <- t(feat.crc.zeller1)
meta.crc.zeller1 = read.table("metadata_sex_age_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE, check.names = F)
#feat.crc.zeller1 = feat.crc.zeller %>% apply(2, function(x) x / 100)
rownames(meta.crc.zeller1)
rownames(meta.crc.zeller1) <-meta.crc.zeller1$Individual
rownames(meta.crc.zeller1)
#data("feat_crc_zeller", package="SIAMCAT")
#data("meta_crc_zeller", package="SIAMCAT")

feat.crc.zeller1[1:3, 1:3]
dim(feat.crc.zeller1)
head(meta.crc.zeller1)

#设置标签
label.crc.zeller1 <- create.label(meta=meta.crc.zeller1, label='study_condition', case='Obesity')
#构建对象
sc.obj1 <- siamcat(feat=feat.crc.zeller1, label=label.crc.zeller1, meta=meta.crc.zeller1)

show(sc.obj1)
#过滤数据集
sc.obj1 <- filter.features(sc.obj1, filter.method = 'abundance', cutoff = 0.00001)
show(sc.obj1)

###################################### Association Testing ################################
sc.obj1 <- check.associations(sc.obj1, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj1, sort.by = 'fc', 
                 panels = c('fc', 'prevalence', 'auroc'))

##################################### Confounder Testing ###################################
check.confounders(sc.obj1, fn.plot = 'kos_confounder_plots.pdf',
                  meta.in = NULL, feature.type = 'filtered')


##################################### Model Building #######################################
#Data Normalization
sc.obj1 <- normalize.features(sc.obj1, norm.method = "log.unit",
                              norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

#Prepare Cross-Validation
sc.obj1 <-  create.data.split(sc.obj1, num.folds = 10, num.resample = 10)

#Model Training
sc.obj1 <- train.model(sc.obj1, method = "lasso")

# get information about the model type
model_type(sc.obj1)
# access the models
models1 <- models(sc.obj1)
models1[[1]]$model

#################################### Make Predictions #######################################
sc.obj1 <- make.predictions(sc.obj1)
pred_matrix1 <- pred_matrix(sc.obj1)

head(pred_matrix1)

################################## Model Evaluation and Interpretation #######################
sc.obj1 <-  evaluate.predictions(sc.obj1)
model.evaluation.plot('vOTU'=sc.obj, 'KOs'=sc.obj1)
model.interpretation.plot(sc.obj1, fn.plot = 'kegg_interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')




################################################################################################
############################ Bacteria abundance table ##############################################

##################################### 加载数据，构架SIAMCAT所需的对象 #####################

feat.crc.zeller2 = read.table("total-7_abundance_table_species_100.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = F)
#feat.crc.zeller2 <- t(feat.crc.zeller2)
meta.crc.zeller2 = read.table("metadata_sex_age_0.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE, check.names = F)
#feat.crc.zeller1 = feat.crc.zeller %>% apply(2, function(x) x / 100)

#data("feat_crc_zeller", package="SIAMCAT")
#data("meta_crc_zeller", package="SIAMCAT")
rownames(meta.crc.zeller2)
rownames(meta.crc.zeller2) <-meta.crc.zeller2$Individual
rownames(meta.crc.zeller2)



feat.crc.zeller2[1:3, 1:3]
dim(feat.crc.zeller2)
head(meta.crc.zeller2)

#设置标签
label.crc.zeller2 <- create.label(meta=meta.crc.zeller2, label='study_condition', case='Obesity')
#构建对象
sc.obj2 <- siamcat(feat=feat.crc.zeller2, label=label.crc.zeller2, meta=meta.crc.zeller2)


show(sc.obj2)
#过滤数据集
sc.obj2 <- filter.features(sc.obj2, filter.method = 'abundance', cutoff = 0.00001)
#sc.obj2 <- filter.features(sc.obj2, filter.method = 'abundance', cutoff = 0.0001)
show(sc.obj2)

###################################### Association Testing ################################
sc.obj2 <- check.associations(sc.obj2, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj2, sort.by = 'fc', 
                 panels = c('fc', 'prevalence', 'auroc'))

##################################### Confounder Testing ###################################
check.confounders(sc.obj2, fn.plot = 'confounder_plots.pdf',
                  meta.in = NULL, feature.type = 'filtered')


##################################### Model Building #######################################
#Data Normalization
sc.obj2 <- normalize.features(sc.obj2, norm.method = "log.unit",
                              norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

#Prepare Cross-Validation
sc.obj2 <-  create.data.split(sc.obj2, num.folds = 10, num.resample = 10)

#Model Training
sc.obj2 <- train.model(sc.obj2, method = "lasso")

# get information about the model type
model_type(sc.obj2)
# access the models
models2 <- models(sc.obj2)
models2[[1]]$model

#################################### Make Predictions #######################################
sc.obj2 <- make.predictions(sc.obj2)
pred_matrix2 <- pred_matrix(sc.obj2)

head(pred_matrix2)

################################## Model Evaluation and Interpretation #######################
sc.obj2 <-  evaluate.predictions(sc.obj2)
model.evaluation.plot('vOTU'=sc.obj, 'KOs'=sc.obj1, 'Bacteria'=sc.obj2)
model.interpretation.plot(sc.obj2, fn.plot = 'kegg_kos_interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

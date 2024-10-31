#对PTV的携带者与非携带者的患病情况，进行逻辑回归
#人群选基因判别的白人，去除基因型缺失大于0.02，sex chromosome aneuploidy，outlier的个体
#协变量选择性别、年龄、年龄平方、测序队列batch和PCA的top10主成分
#蛋白以ABCDE举例，可以替换成你的蛋白

library(data.table)
library(tidyverse)

#协变量数据
covariate1 <- fread("/share/home/zhangjunyu/Project/240925_HF/Data/Covariates/covariate_before_inputation.csv")
covariate2 <- fread("/share/home/zhangjunyu/Rawdata/UKB_GWAS/GWAS_cov.csv")
##人群QC
covariate2 <- covariate2 %>% filter(Genetic_ethnic_grouping == "Caucasian")
covariate2 <- covariate2 %>% filter(Sex_chromosome_aneuploidy != "Yes")
covariate2 <- covariate2 %>% filter(Outliers_for_heterozygosity_or_missing_rate != "Yes")
covariate2 <- covariate2 %>% filter(Missingness <= 0.02)
##协变量提取及合并
covariate1 <- covariate1 %>% select(eid, age)
covariate1 <- covariate1 %>% mutate(age_square = age^2)
covariate2 <- covariate2 %>% mutate(genotyping_array = ifelse(grepl("Batch_", Genotype_measurement_batch), 1, 0))
covariate2 <- covariate2 %>% select(eid, genotyping_array, Genetic_sex, Genetic_principal_components_Array_1, Genetic_principal_components_Array_2, 
Genetic_principal_components_Array_3, Genetic_principal_components_Array_4, Genetic_principal_components_Array_5, Genetic_principal_components_Array_6,
Genetic_principal_components_Array_7, Genetic_principal_components_Array_8, Genetic_principal_components_Array_9, Genetic_principal_components_Array_10)
covariate <- inner_join(covariate1, covariate2, by = "eid")
covariate <- covariate %>% mutate(age_square = as.integer(age_square), genotyping_array = factor(genotyping_array), Genetic_sex = factor(Genetic_sex))

#IV数据
setwd("/share/home/linmiao/bio/2023_AS_EC/LIMA1/UKB_PTVs/ABCDE")
IV_status <- fread("ABCDE_LOF_variations.raw")
##IV加工
IV_status <- IV_status[,-2:-6]
IV_status <- IV_status %>% rename(eid = FID)
names(IV_status) <- gsub("(.*)_.*", "\\1",names(IV_status))

#第一种结局：慢性肾病
##读取结局数据
phenotype <- fread("/share/home/linmiao/bio/2023_AS_EC/LIMA1/UKB_PTVs/Disease_outcomes/ASCVD_outcomes.csv")
##LOFTEE
LOFTEE_list <- fread("ABCDE_IVs_LOFTEE.txt", header = F)
variation_list <- unique(unlist(c("eid", LOFTEE_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)

##AlphaMissense
AlphaMissense_list <- fread("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/genetic_intervention/ABCDE/ABCDE_IVs_AlphaMissense.txt", header = F)
variation_list <- unique(unlist(c("eid", AlphaMissense_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)


##SIFT
SIFT_list <- fread("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/genetic_intervention/ABCDE/ABCDE_IVs_SIFT.txt", header = F)
variation_list <- unique(unlist(c("eid", SIFT_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)

##LOFTEE+AlphaMissense
variation_list <- unique(unlist(c("eid", LOFTEE_list, AlphaMissense_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)

##LOFTEE+SIFT
variation_list <- unique(unlist(c("eid", LOFTEE_list, SIFT_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)

##AlphaMissense+SIFT
variation_list <- unique(unlist(c("eid", SIFT_list, AlphaMissense_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)

##LOFTEE+AlphaMissense+SIFT
variation_list <- unique(unlist(c("eid", LOFTEE_list, AlphaMissense_list, SIFT_list)))
variation_status <- IV_status %>% select(any_of(variation_list))
variation_status$sum <- rowSums(variation_status[, -1], na.rm =T)  # -1表示排除第一列 FID
variation_status <- variation_status %>% mutate(PAVcarrier = ifelse(sum == 0, 0, 1))
variation_status <- variation_status %>% select(eid, PAVcarrier)
variation_status <- left_join(variation_status, phenotype, by = "eid")
variation_status <- variation_status %>% mutate(target_y = ifelse(is.na(target_y), 1, target_y))
variation_status <- variation_status %>% rename(outcome = "target_y")
variation_status <- variation_status %>% mutate(PAVcarrier = factor(PAVcarrier), outcome = factor(outcome))
data_all <- inner_join(variation_status, covariate, by = "eid")
model <- glm(outcome ~ PAVcarrier + Genetic_sex + age + age_square + Genetic_principal_components_Array_1 + Genetic_principal_components_Array_2 +
                       Genetic_principal_components_Array_3 + Genetic_principal_components_Array_4 + Genetic_principal_components_Array_5 + 
                       Genetic_principal_components_Array_6 + Genetic_principal_components_Array_7 + Genetic_principal_components_Array_8 +
                       Genetic_principal_components_Array_9 + Genetic_principal_components_Array_10, data = data_all, family = binomial)
summary(model)
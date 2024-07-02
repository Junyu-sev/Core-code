#基于最大化Youden Index获取每个蛋白的cutoff数值，和高低浓度组未来患疾病的HR和P-val

#加载R包
library(tidyverse)
library(data.table)
library(survival)
library(mice)


#路径定义-----------------------------------------------------------------------
proteomics_path <- "/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/Proteomics/"                                                           #含有蛋白组学数据的路径
covariate_path <- "/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/covariate/"                                                            #含有协变量数据的路径
disease_def_path <- "/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/DKD_outcome/DKD_outcome_complication_baseline_removed.csv"           #在这里更改疾病定义的路径
cox_path <- "/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/Cox/240701_Cox/DKD"
cutoff_path <- '/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/KM/240701_cutoff/DKD_cutoff'                                          #在这里更改cutoff输出的路径


##############################################################################
###                                                                        ###
###                                KM生存分析                               ###
###                                                                        ###
##############################################################################

print("开始：KM生存分析")

#通过计算youden指数，获取每个蛋白的cutoff数值
target_df <- fread(disease_def_path, select = c('eid', 'target_y', 'BL2Target_yrs'))
pro_df <- fread (paste(proteomics_path, '231008_data_olink_instance_0.csv',sep=""))
cox_df <- fread(paste(cox_path, ".csv",sep=""))
pro_list <- cox_df %>% filter(p_val_bfi < 0.05)
pro_list <- pro_list$Pro_code

mydf <- inner_join(target_df, pro_df, by = "eid")

cox_df <- cox_df %>% mutate(Factor = ifelse(Coef < 0, "Protective", "Risk"))

get_Youden <- function(mydf, cutoff, Factor){
  if(Factor == "Risk"){
    mydf <- mydf %>% mutate(pred_y = ifelse(target_pro > cutoff, 1, 0))
  }else{
    mydf <- mydf %>% mutate(pred_y = ifelse(target_pro > cutoff, 0, 1))
  }
  TP <- sum(mydf$target_y == 1 & mydf$pred_y == 1)
  FP <- sum(mydf$target_y == 0 & mydf$pred_y == 1)
  FN <- sum(mydf$target_y == 1 & mydf$pred_y == 0)
  TN <- sum(mydf$target_y == 0 & mydf$pred_y == 0)
  Sensitivity <- TP/(TP + FN)
  Specificity <- TN/(TN + FP)
  Youden_index <- Sensitivity + Specificity - 1
  return(Youden_index)
}

cutoff_result <- c()

i=0
for(f in pro_list){
 i=i+1
 print(i)
 tmp_df <- mydf %>% select(any_of(c("target_y", f)))
 names(tmp_df)[2] <- "target_pro"
 tmp_df <- na.omit(tmp_df)
 bin <- (max(tmp_df$target_pro, na.rm = T) - min(tmp_df$target_pro, na.rm = T))/9999
 cutoff_lst <- seq(min(tmp_df$target_pro, na.rm = T), max(tmp_df$target_pro, na.rm = T), by = as.numeric(bin))
 Youden <- tibble(cutoff_lst)
 Factor <- cox_df %>% filter(Pro_code == f) %>% {.$Factor}
 Youden$Youden_index <- sapply(cutoff_lst, function(cutoff) get_Youden(tmp_df, cutoff, Factor))
 cutoff_add <- Youden[which(Youden$Youden_index == max(Youden$Youden_index, na.rm = T)),1]
 cutoff_result <- c(cutoff_result, as.numeric(cutoff_add[1,1]))
}

result <- tibble(pro_list, cutoff_result)
names(result)[1] <- "Pro_code"
write.csv(result, paste(cutoff_path, '.csv', sep = ""), row.names = F)

#根据每个蛋白的cutoff值，计算相应的HR和p-val
load(file = paste(covariate_path,"DKD_imp_5.RData",sep=""))
imp_cov_data <- complete(imp.cov, "long", include = TRUE)

target_df <- fread(disease_def_path, select = c('eid', 'target_y', 'BL2Target_yrs'))
pro_df <- fread (paste(proteomics_path, '231008_data_olink_instance_0.csv',sep=""))
mydf <- inner_join(target_df, pro_df, by = 'eid')
cutoff <- fread(paste(cutoff_path, ".csv", sep = ""))

imp_cov_data$BL2Target_yrs <- rep(mydf$BL2Target_yrs, imp.cov$m + 1)
imp_cov_data$target_y <- rep(mydf$target_y, imp.cov$m + 1)
diet_cols <- names(imp_cov_data)[grep("_score", names(imp_cov_data))]
imp_cov_data[, diet_cols] <- lapply(imp_cov_data[, diet_cols], function(f) as.integer(as.character(f)))
imp_cov_data$Diet <- rowSums(imp_cov_data[, diet_cols])
imp_cov_data <- imp_cov_data[, !names(imp_cov_data) %in% diet_cols]

myout_df <- tibble()
i = 0
for (pro_f in cutoff$Pro_code){
  i = i + 1
  print(i)
  imp_cov <- imp_cov_data
  imp_cov$target_pro <- rep(mydf[,get(pro_f)], imp.cov$m + 1)
  imp_cov <- imp_cov %>% filter(!is.na(target_pro))
  imp_cov <- imp_cov %>% mutate(target_pro = ifelse(target_pro > cutoff$cutoff_result[which(cutoff$Pro_code == pro_f)], 1, 0))
  imp_cov$target_pro <- factor(imp_cov$target_pro)
  imp_cov$.id <- rep(1:(nrow(imp_cov)/6), imp.cov$m + 1)
  imp.cov.new <- as.mids(imp_cov)
  tryCatch({
    fit <- with(imp.cov.new, coxph(Surv(BL2Target_yrs, target_y) ~ target_pro + age + sex + Ethnic + Qualification + Townsend + income + Smoking + Alcohol + IPAQ  + BMI + HbA1c + hypertension + Diabets_time + Treatment + Diet))
    Coef <- summary(pool(fit))$estimate[1]
    HR <- exp(Coef)
    sd.coef <- summary(pool(fit))$std.error[1]
    pval<- summary(pool(fit))$p.value[1]
    myout = data.frame(pro_f, Coef, sd.coef, HR, pval)
    myout_df = bind_rows(myout_df,myout)
  }, error = function(e){
    print(pro_f)
  })
  }
names(myout_df) <- c('Pro_code', 'Coef', 'sd.Coef', 'HR', 'HR_p_val')
myout_df$p_val_bfi <- p.adjust(myout_df$HR_p_val, method ="bonferroni")
myout_df <- myout_df %>% rowwise() %>% mutate(HR_ci = paste(round(HR, 2), " [", round(exp(Coef - 1.96*sd.Coef), 2), "-", round(exp(Coef + 1.96*sd.Coef),2), "]", sep = ""))
myout_df <- full_join(cutoff, myout_df, by = "Pro_code")
myout_df$p_val_bfi <- signif(myout_df$p_val_bfi,3)
write.csv(myout_df, paste(cutoff_path, '_HR.csv', sep = ""), row.names = F)

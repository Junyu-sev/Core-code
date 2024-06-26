#多重插补的教程参考https://www.bookdown.org/rwnahhas/RMPH/mi-fitting.html
#多重插补的教程参考https://zhuanlan.zhihu.com/p/584937019
#多重插补的教程参考https://stefvanbuuren.name/fimd/
#多重插补的教程参考https://zhuanlan.zhihu.com/p/584937019

library(data.table)
library(tidyverse)
library(survival)
library(mice)
#library(miceadds)
#library(howManyImputations)
#source("Functions_rmph.R")

dpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Data/'
outpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/cox_with_covariate_imputed/'

cov_df <- fread(paste(dpath,"Raw_data/covariate_before_inputation.csv",sep=""))
cov_df <- as_tibble(cov_df)

char_columns <- as.logical(sapply(cov_df, is.character))
cov_df[, char_columns] <- lapply(cov_df[, char_columns], as.factor)
num_columns <- c("sex","Cholesterol_lowering_medication","Blood_pressure_medication","hypertension","Illnesses_of_family_stroke_and_heart_disease")
cov_df[, num_columns] <- lapply(cov_df[, num_columns], as.factor)

#转为有序因子变量
cov_df$income <- ordered(cov_df$income, 
      levels = c("Less than 18 000","18 000 to 30 999","31 000 to 51 999","52 000 to 100 000","Greater than 100 000"))
cov_df$IPAQ <- ordered(cov_df$IPAQ, 
                            levels = c("low","moderate","high"))
cov_df$Overall_health_rating <- ordered(cov_df$Overall_health_rating, 
                            levels = c("Poor","Good","Fair","Excellent"))                 
m1_f_lst <- names(cov_df)[-1]

#多重插补
imp.cov <- mice(cov_df[,-1],
                   seed  = 2024,
                   m     = 5,
                   print = F)
save.image(file = paste(outpath,"imp.RData",sep=""))


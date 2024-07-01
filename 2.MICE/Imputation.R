#多重插补的教程参考https://www.bookdown.org/rwnahhas/RMPH/mi-fitting.html
#多重插补的教程参考https://zhuanlan.zhihu.com/p/584937019
#多重插补的教程参考https://stefvanbuuren.name/fimd/

#R包加载
library(data.table)
library(tidyverse)
library(mice)

#路径设置
dpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Data/'                                   #协变量读取路径
outpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/cox_with_covariate_imputed/'    #插补结果输出路径，将以.RData格式输出

cov_df <- fread(paste(dpath,"covariate_before_inputation.csv",sep=""))
cov_df <- as_tibble(cov_df)

#在插补前，把字符型变量变为因子变量；这里有些变量存储形式是1（表示有）/0（表示没有），所以有些形式为数值的列也要转为因子变量
char_columns <- as.logical(sapply(cov_df, is.character))
cov_df[, char_columns] <- lapply(cov_df[, char_columns], as.factor)
#哪些列的存储形式是1（表示有）/0（表示没有），例如Cholesterol_lowering_medication中“1”表示有降脂治疗，“0”表示没有
num_columns <- c("sex","Cholesterol_lowering_medication","Blood_pressure_medication","hypertension","Illnesses_of_family_stroke_and_heart_disease") 
cov_df[, num_columns] <- lapply(cov_df[, num_columns], as.factor)

#部分因子变量转为有序因子变量（这些变量的取值是有顺序，有高低之分的）
cov_df$income <- ordered(cov_df$income, 
      levels = c("Less than 18 000","18 000 to 30 999","31 000 to 51 999","52 000 to 100 000","Greater than 100 000"))
cov_df$IPAQ <- ordered(cov_df$IPAQ, 
                            levels = c("low","moderate","high"))
cov_df$Overall_health_rating <- ordered(cov_df$Overall_health_rating, 
                            levels = c("Poor","Good","Fair","Excellent"))                 
# 去掉eid列
m1_f_lst <- names(cov_df)[-1]

#多重插补，插补5次
imp.cov <- mice(cov_df[,-1],
                   seed  = 2024,
                   m     = 5,
                   print = F)
save.image(file = paste(outpath,"imp.RData",sep=""))


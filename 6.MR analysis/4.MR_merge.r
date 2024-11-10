setwd("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/MR/UKB_Olink_Finngen_DKD_1E-5/UKB_Olink_Finngen_Diabetic_nephropathy")

library(data.table)
library(tidyverse)

files <- list.files()
MR <- tibble()

for(file in files){
    dat <- fread(file)
    if(nrow(dat)!=1){
        MR <- bind_rows(MR, dat[3,])
    }else{
        MR <- bind_rows(MR, dat)
    }
}
MR <- MR[,-1:-2]
MR <- MR[order(MR$pval),]

write.csv(MR, "/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/MR/UKB_Olink_Finngen_DKD_1E-5/CSV_integration/UKB_Olink_Finngen_Diabetic_nephropathy.csv", row.names = F)
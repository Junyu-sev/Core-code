### 需要满足以下要求才会纳入：
### 1.有蛋白组学数据
### 2.基线前不患病

#加载R包
library(tidyverse)
library(data.table)

#疾病定义-----------------------------------------------------------------------
target_disease <- "Stroke"                                   #在这里更改疾病名称
disease_ICD10 <- "I60|I61|I62|I63|I64|I65|I66|I67|I68|I69"   #在这里更改疾病的ICD10编码，用“或者”的符号“|”进行连接
disease_OPSC4 <- ""                                          #在这里更改疾病手术的OPSC编码，用“或者”的符号“|”进行连接

#路径定义-----------------------------------------------------------------------
proteomics_path <- "/share/home/zhangjunyu/Project/Proteomic_analysis/Data/Proteomics/"                                       #含有蛋白组学数据的路径
diagnose_path <- "/share/home/zhangjunyu/Project/Proteomic_analysis/Data/Diagnosis_info/"                                     #含有疾病诊断信息的路径
disease_def_path <- "/share/home/zhangjunyu/Project/Proteomic_analysis/Data/Disease_outcomes/Stroke/"                         #在这里更改疾病定义输出的路径

##############################################################################
###                                                                        ###
###                  疾病定义 纳入排除 获取生存分析数据                       ###
###                                                                        ###
##############################################################################

print("开始：疾病定义 纳入排除 获取生存分析数据")

#获取各种数据-------------------------------------------------------------------
##获取蛋白组学数据的人的索引，后续数据排除没有蛋白组学数据的个体
proteomics_eid <- fread(paste(proteomics_path, "231008_data_olink_instance_0.csv", sep = ""), select = "eid")
proteomics_eid <- unlist(proteomics_eid)
##获取出生年份，当自我报告用年龄表示疾病时间时，能转为完整时间
date_birth <- fread(paste(diagnose_path, "year of birth.csv", sep = ""))
date_birth <- subset(date_birth, date_birth$eid %in% proteomics_eid)
##获取入组时间
attending_time <- fread(paste(covariate_path, "covariate.csv", sep = ""),select = c("eid", "attending_time"))
attending_time<- subset(attending_time, eid %in% proteomics_eid)
attending_time$attending_time <- ymd(attending_time$attending_time)
##获取住院数据
inpatient <- fread(paste(diagnose_path, "inpatient diagnose (date included).csv", sep = ""))
inpatient <- subset(inpatient, inpatient$Participant.ID %in% proteomics_eid)
inpatient <- rename(inpatient,eid = Participant.ID)
##获取自我报告数据
self_reported <- fread(paste(diagnose_path, "self_reported (date included).csv", sep = ""))
self_reported <- subset(self_reported, self_reported$Participant.ID %in% proteomics_eid)
self_reported <- rename(self_reported, eid = Participant.ID)
##获取死亡记录数据
death <- fread(paste(diagnose_path, "death disease (date included).csv", sep = ""))
death <- subset(death, death$Participant.ID %in% proteomics_eid)
death <- rename(death, eid = Participant.ID)
##获取primary care数据
primary_care <- fread(paste(diagnose_path, "primary care (date included).csv", sep = ""), 
     select = c("Participant ID","Date clinical code was entered","Read v2 code","CTV3 (Read v3) code"))
primary_care  <- subset(primary_care , primary_care $`Participant ID` %in% proteomics_eid)
primary_care <- rename(primary_care, eid = `Participant ID`)
##获取手术事件
if(disease_OPSC4 != ""){
  operation <- fread(paste(diagnose_path, "operation_record (date included).csv", sep = ""))
  operation <- subset(operation, operation$`Participant ID` %in% proteomics_eid)
  operation <- rename(operation, eid = `Participant ID`)
}

#住院诊断数据的处理-------------------------------------------------------------
##将疾病编码替换为单纯的ICD10编码形式
inpatient <- inpatient[grep(pattern = disease_ICD10, inpatient$Diagnoses...ICD10),]
diagnose <- strsplit(unlist(inpatient[,2]),split = "\\|")
diagnose <- lapply(diagnose, function(x) substr(x, start = 1, stop = 5))
diagnose <- lapply(diagnose, function(x) gsub(pattern = " .","",x))
diagnose <- sapply(diagnose, function(x) paste(x,collapse = "|"))
inpatient$Diagnoses...ICD10 <- diagnose

##提取想要的疾病患者
get_disease <- function(x){
  disease <- strsplit(x[2], split = "\\|")
  disease <- unlist(disease)
  index <- grep(pattern = disease_ICD10, disease) + 2
  x[3] <- paste(x[index], collapse = "|")
  disease <- disease[grep(pattern = disease_ICD10, disease)]
  disease <- paste(disease, collapse = "|")
  x[2] <- disease
  x <- x[c(1, 2, 3)]
  return(x)
}
inpatient <- apply(inpatient, 1, get_disease)
inpatient <- as_tibble(t(inpatient))
names(inpatient) <- c("eid", "inpatient", "date_inpatient")

#self_reported数据的处理--------------------------------------------------------
##选取有编码转化关系的non-cancer疾病（因此，cancer疾病的患者统计是有问题的，需要调整代码）
cols <- rep("Non.cancer.illness.year.age.first.occurred...Instance.0...Array.",34)
cols <- paste(cols,0:33,sep = "")
self_reported <- self_reported %>% select("eid","Non.cancer.illness.code..self.reported...Instance.0",any_of(cols))
disease <- strsplit(unlist(self_reported[,2]),split = "\\|")

##获取编码转化的方式，用coding保存
coding_self_reported <- fread(paste(diagnose_path, "ICD10mapping/coding_self_reported.tsv", sep = ""))
coding_self_reported <- coding_self_reported %>% filter(coding > 0)
coding_self_reported <- coding_self_reported[,1:2]
coding_self_reported_to_ICD10 <- fread(paste(diagnose_path, "ICD10mapping/coding_self_reported_to_ICD10.tsv", sep = ""))
coding <- merge(coding_self_reported,coding_self_reported_to_ICD10,by="coding")
coding <- coding[,2:3]
names(coding) <- c("Disease","ICD10")
coding[,1] <- gsub("\\(","\\\\(",coding$Disease)#确保R的正则表达式是ok的
coding[,1] <- gsub("\\)","\\\\)",coding$Disease)

##将自我报告的编码替换为ICD10编码
replace_self_reported_with_ICD10 <- function(disease) {
  for (i in 1:nrow(coding)) {
    disease <- gsub(
      paste0("^", coding$Disease[i], "$"), coding$ICD10[i], disease)}
  return(disease)
}#定义了自我报告向ICD10转化的函数，下面实现转化
disease <- lapply(disease, replace_self_reported_with_ICD10)
disease <- sapply(disease, function(x) paste(x, collapse = "|"))
self_reported[,2] <- disease

##挑选包含疾病的患者，认为是基线前患病的
self_reported <- self_reported[grep(pattern = disease_ICD10, self_reported$Non.cancer.illness.code..self.reported...Instance.0),]


#death数据的处理----------------------------------------------------------------
##death诊断数据仅保留ICD10编码，合并成一列
disease <- apply(death[,2:32], 2, function(col) substr(col, start = 1, stop = 5))
disease <- apply(disease, 2, function(col) gsub(pattern = " .","",col))
disease <- apply(disease,1,function(row) unique(row))
disease <- sapply(disease, function(x) paste(x,collapse = "|"))
disease <- gsub(pattern = "NA\\|","",disease)
disease <- gsub(pattern = "\\|\\|","\\|",disease)
disease <- gsub(pattern = "NA","",disease)
disease <- gsub(pattern = "\\|$","",disease)
death[,2] <- disease
death <- death[,c(1,2,34)]

##获取目标疾病诊断的个体
names(death) <- c("eid","death","death_date")
death <- death[grep(pattern = disease_ICD10, death),]

##death数据有些个体数据存在重复（原始数据中就有重复），处理一下
death <- split(death, death$eid)
death <- lapply(death,function(x){
  x[1,2] <- paste(unlist(x[,2]),collapse = "|")
  x <- x[1,]
  y <- strsplit(as.character(x[1,2]), split = "\\|")
  y <- unlist(y)
  y <- unique(y)
  x[1,2] <- paste(y,collapse = "|")
  return(x)
})
death <- do.call(rbind, death)
get_disease_death <- function(x){
  disease <- strsplit(x[2], split = "\\|")
  disease <- unlist(disease)
  disease <- disease[grep(pattern = disease_ICD10, disease)]
  disease <- paste(disease, collapse = "|")
  x[2] <- disease
  return(x)
}
death <- apply(death, 1, get_disease_death)
death <- as_tibble(t(death))
names(death) <- c("eid", "death", "date_death")

#primary care数据的处理---------------------------------------------------------
##获取primary care疾病名称和ICD10对应关系，进行map
coding_read2_to_ICD10 <- fread(paste(diagnose_path, "ICD10mapping/coding_Read2_to_ICD10.tsv", sep = ""))
coding_read3_to_ICD10 <- fread(paste(diagnose_path, "ICD10mapping/coding_Read3_to_ICD10.tsv", sep = ""))
primary_care <- merge(primary_care,coding_read2_to_ICD10,by.x = "Read v2 code", by.y = "coding",all.x = T, all.y = F)
primary_care <- merge(primary_care,coding_read3_to_ICD10,by.x = "CTV3 (Read v3) code", by.y = "coding",all.x = T, all.y = F)
primary_care <- primary_care[,-1:-2]
names(primary_care) <- c("eid","Date","ICD10_2","ICD10_3")

##获取目标疾病的患者
primary_care <- primary_care[c(grep(pattern = disease_ICD10, primary_care$ICD10_2), grep(pattern = disease_ICD10, primary_care$ICD10_3)),]
primary_care <- primary_care[order(primary_care$eid),]
primary_care <- subset(primary_care, primary_care$Date != "")

#因为primary care数据每个个体可能不只有一行数据，所以要把有多行数据的个体处理一下
primary_care <- split(primary_care, primary_care$eid)

format_change_primary_care <- function(x){
  Participant.ID <- x[1,1]
  Date <- paste(unlist(x[,2]), collapse = "|")
  if (!is.na(x[1,3])){
    ICD10 <- paste(unlist(x[,3]), collapse = "|")
    x <- data.frame(Participant.ID, ICD10, Date)
    return(x)
  }
  if (!is.na(x[1,4])){
    ICD10 <- paste(unlist(x[,4]), collapse = "|")
    x <- data.frame(Participant.ID, ICD10, Date)
    return(x)
  }
}
primary_care <- lapply(primary_care, format_change_primary_care)
primary_care <- do.call(rbind, primary_care)
names(primary_care) <- c("eid", "primariy_care", "date_primary_care")

#手术事件的处理-----------------------------------------------------------------
if(disease_OPSC4 != ""){
  operation <- operation[grep(pattern = disease_OPSC4, operation$`Operative procedures - OPCS4`),]
  diagnose <- strsplit(unlist(operation[,2]),split = "\\|")
  diagnose <- lapply(diagnose, function(x) substr(x, start = 1, stop = 5))
  diagnose <- lapply(diagnose, function(x) gsub(pattern = " .","",x))
  diagnose <- sapply(diagnose, function(x) paste(x,collapse = "|"))
  operation$`Operative procedures - OPCS4` <- diagnose
  rm(diagnose)
  get_disease_operation <- function(x){
    disease <- strsplit(x[2], split = "\\|")
    disease <- unlist(disease)
    index <- grep(pattern = disease_OPSC4, disease) + 2
    x[3] <- paste(x[index], collapse = "|")
    disease <- disease[grep(pattern = disease_OPSC4, disease)]
    disease <- paste(disease, collapse = "|")
    x[2] <- disease
    x <- x[c(1, 2, 3)]
    return(x)
}
  operation <- apply(operation, 1, get_disease_operation)
  operation <- as_tibble(t(operation))
  names(operation) <- c("eid", "operation", "date_operation")
}

#合并数据,排除自我报告这种疾病的个体（认为自我报告均在基线前患病）-------------------
diagnose_date <- merge(inpatient[,c(1,3)], death[,c(1,3)], by = "eid", all = T)
diagnose_date <- merge(diagnose_date, primary_care[,c(1,3)], by = "eid", all = T)
diagnose_date <- subset(diagnose_date, !diagnose_date$eid %in% self_reported$eid)

if(disease_OPSC4 != ""){
  diagnose_date <- merge(operation[,c(1,3)], inpatient[,c(1,3)], by = "eid", all = T)
  diagnose_date <- merge(diagnose_date, death[,c(1,3)], by = "eid", all = T)
  diagnose_date$Participant.ID <- as.integer(diagnose_date$eid)
  diagnose_date <- merge(diagnose_date, primary_care[,c(1,3)], by = "eid", all = T)
  diagnose_date <- subset(diagnose_date, !diagnose_date$eid %in% self_reported$eid)
}


diagnose_date <- apply(diagnose_date, 1, function(x){
  date <- x[-1]
  date <- paste(date, collapse = "|")
  x[2] <- date
  x <- x[1:2]
  return(x)
})
diagnose_date <- t(diagnose_date)
diagnose_date <- as_tibble(diagnose_date)
names(diagnose_date) <- c("eid", "diagnose_date")

diagnose_date$diagnose_date <- gsub("\\|NA", "", diagnose_date$diagnose_date)
diagnose_date$diagnose_date <- gsub("NA\\|", "", diagnose_date$diagnose_date)

max_cols <- max(sapply(strsplit(diagnose_date$diagnose_date, "\\|"), length))
diagnose_date <- separate(diagnose_date, diagnose_date, into = paste0("date_", 1:max_cols), sep = "\\|")
date_cols <- names(diagnose_date)[2:ncol(diagnose_date)]
diagnose_date <- mutate(diagnose_date,across(all_of(date_cols), as.Date))
diagnose_date <- merge(diagnose_date, attending_time, by = "eid")

#计算基线到诊断的时间跨度,提供疾病组数据
diagnose_date1 <- diagnose_date[,c(-1, -ncol(diagnose_date))]
diagnose_date1 <- apply(diagnose_date1, 2, function(x){
  x <- as.Date(x)
  x <- tibble(x, as.Date(diagnose_date$attending_time))
  x <- x[,1] - x[,2]
  x <- as.numeric(unlist(x))/365
  return(x)
})
diagnose_date1 <- apply(diagnose_date1, 1, min, na.rm = T)#获取最早出现事件的时间
diagnose_date[2] <- diagnose_date1
diagnose_date <- diagnose_date[,1:2]
baseline_diagnose_date <- diagnose_date %>% filter(date_1 <= 0)
names(baseline_diagnose_date) <- c("eid", "BL2Target_yrs")
diagnose_date <- diagnose_date %>% filter(date_1 > 0)
diagnose_date <- diagnose_date %>% filter(!is.infinite(date_1))
diagnose_date$target_y <- 1
names(diagnose_date) <- c("eid", "BL2Target_yrs", "target_y")

#制作对照数据，去掉自我报告，去掉基线前诊断，去掉基线后诊断
control <- attending_time
control <- subset(control, !control$eid %in% self_reported$eid)
control <- subset(control, !control$eid %in% baseline_diagnose_date$eid)
control <- subset(control, !control$eid %in% diagnose_date$eid)

control <- as_tibble(control)
##获取死亡时间
deathdate <- fread(paste(diagnose_path, "death disease (date included).csv", sep = ""))
deathdate <- deathdate[,c(1,34)]
names(deathdate) <- c("eid", "date_death")
deathdate <- distinct(deathdate)

#计算时间跨度
control <- merge(control, deathdate, by = "eid", all.x = T)
control$record_end <- "2022-10-31"#这是住院、死亡记录数据最后更新的时间点
control <- as_tibble(control)
control$attending_time <- as.Date(control$attending_time)
control$date_death <- as.Date(control$date_death)
control$record_end <- as.Date(control$record_end)
control$date_death <- as.numeric(control$date_death - control$attending_time)/365
control$record_end <- as.numeric(control$record_end - control$attending_time)/365
control <- control %>% mutate(BL2Target_yrs = ifelse(!is.na(date_death), date_death, record_end))#很奇怪的一件事是，有些个体的死亡记录
control$target_y <- 0
control <- control %>% select(c("eid", "BL2Target_yrs", "target_y"))

Disease_outcomes <- rbind(diagnose_date, control)
Disease_outcomes <- Disease_outcomes[order(Disease_outcomes$eid),]
write.csv(Disease_outcomes, paste(disease_def_path , target_disease, "_outcomes.csv", sep = ""), row.names = F)
save.image(file = paste(disease_def_path , target_disease, ".RData", sep = ""))

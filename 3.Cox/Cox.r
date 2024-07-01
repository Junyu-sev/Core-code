#对每种UKB蛋白，分别做疾病的cox比例风险模型分析，获取每个蛋白的HR和校正p-val
#cox协变量要根据具体疾病进行调整，下面采用的协变量应用于IHD（缺血性心脏病）

#加载R包
library(tidyverse)
library(data.table)
library(survival)
library(mice)
library(ggplot2)
library(ggrepel)


#路径定义-----------------------------------------------------------------------
proteomics_path <- "/share/home/zhangjunyu/Project/IHD_Proteomic_analysis/Data/Proteomics/"                       #含有蛋白组学数据的路径
covariate_path <- "/share/home/zhangjunyu/Project/IHD_Proteomic_analysis/Data/Covariates/"                        #含有协变量数据的路径
disease_def_path <- "/share/home/zhangjunyu/Project/IHD_Proteomic_analysis/Data/Disease_outcomes/IHD/IHD.csv"     #在这里更改疾病定义的路径
cox_path <- "/share/home/zhangjunyu/Project/IHD_Proteomic_analysis/Result/Cox/IHD"                                #在这里更改输出cox结果的路径


##############################################################################
###                                                                        ###
###                        cox分析找寻风险蛋白因素                          ###
###                                                                        ###
##############################################################################

print("开始：cox分析找寻风险蛋白因素")

#加载多重插补的结果
load(file = paste(covariate_path,"imp.RData",sep=""))
imp_cov_data <- complete(imp.cov, "long", include = TRUE)

#合并蛋白组学数据和疾病定义数据
pro_df <- fread (paste(proteomics_path, '231008_data_olink_instance_0.csv',sep=""))
pro_f_lst <- names(pro_df)[-1]
target_df <- fread(disease_def_path, select = c('eid', 'target_y', 'BL2Target_yrs'))
mydf <- inner_join(target_df, pro_df, by = 'eid')

#将合并数据导入多重插补的结果
index <- which(cov_df$eid %in% mydf$eid)
imp_cov_data <- imp_cov_data %>% filter(.id %in% index)
imp_cov_data$BL2Target_yrs <- rep(mydf$BL2Target_yrs, imp.cov$m + 1)
imp_cov_data$target_y <- rep(mydf$target_y, imp.cov$m + 1)

#循环对每个蛋白进行cox分析，获取HR和p-val
myout_df <- tibble()
pro_out_lst <- c()  
i = 0
for (pro_f in pro_f_lst){
  i = i + 1
  print(i)
  imp_cov <- imp_cov_data
  # 在插补数据中加入每个个体当前分析蛋白的浓度
  imp_cov$target_pro <- rep(mydf[,get(pro_f)], imp.cov$m + 1)
  # 删除蛋白浓度缺失的个体
  imp_cov <- imp_cov %>% filter(!is.na(target_pro))
  imp_cov$.id <- rep(1:(nrow(imp_cov)/6), imp.cov$m + 1)
  imp.cov.new <- as.mids(imp_cov)
  tryCatch({
    #协变量要根据具体需求调整
    fit <- with(imp.cov.new, coxph(Surv(BL2Target_yrs, target_y) ~ target_pro + age + sex + Ethnic + Qualification + Townsend + income + assessment_centre + Smoking + Alcohol + IPAQ + Overall_health_rating + BMI + Weight + Height + Waist_circumference + SBP + DBP + Cholesterol + HDL_cholesterol + Triglycerides + HbA1c + hypertension + Diabetes + Cholesterol_lowering_medication + Blood_pressure_medication + Illnesses_of_family_stroke_and_heart_disease))
    Coef <- summary(pool(fit))$estimate[1]
    HR <- exp(Coef)
    sd.coef <- summary(pool(fit))$std.error[1]
    pval<- summary(pool(fit))$p.value[1]
    myout = data.frame(Coef, sd.coef, HR, pval)
    myout_df = bind_rows(myout_df,myout)
    pro_out_lst <- c(pro_out_lst,pro_f) 
  }, error = function(e){
    print(pro_f)
  })
  }
  names(myout_df) <- c('Coef', 'sd.Coef', 'HR', 'HR_p_val')
  myout_df$Pro_code <- pro_out_lst
  myout_df <- myout_df %>% mutate(HR_p_val = ifelse(is.na(HR_p_val),1,HR_p_val)) 
  #进行p值校正
  myout_df$p_val_fdr <- p.adjust(myout_df$HR_p_val, method ="BH")
  myout_df$p_val_bfi <- p.adjust(myout_df$HR_p_val, method ="bonferroni")
  write.csv(myout_df,paste(cox_path,".csv",sep=""), row.names = F)

#对于所有蛋白绘制火山图
myout_df$selectedpro <- ifelse(myout_df$p_val_bfi<0.05, myout_df$Pro_code,NA)
myout_df$selected <- ifelse(myout_df$p_val_bfi<0.05, "Risk protein","Non-risk protein")
p1 <- ggplot(myout_df,aes(HR,-log10(p_val_bfi),
	color = factor(selected),
	size = factor(selected)))+  
      geom_point()+
      labs(x="HR",
           y=expression(-Log[10]*" (p value)"))+
      theme_classic(base_size = 15)+
      scale_color_manual(values = c('grey','darkred'))+
      scale_size_manual(values = c(1,2))+
      geom_text_repel(aes(label=selectedpro), color="black",size=3,
                      box.padding=unit(0.5, "lines"), 
                      point.padding=NA, 
                      segment.colour = "black")+
      theme(legend.title = element_blank(),
      	legend.position = c(0.85,0.95),
      	legend.background = element_rect(fill='transparent'))+
      geom_hline(yintercept = -log10(0.05),linetype=2,cex=1)+
      geom_vline(xintercept = 1,linetype=2,cex=1) 
ggsave(paste(cox_path,".pdf",sep=""),p1)

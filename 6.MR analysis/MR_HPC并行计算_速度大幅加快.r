#用30核来计算所有疾病的MR，5E-8
library(parallel)
library(genetics.binaRies)
library(ieugwasr)
library(gwasglue)
library(gwasvcf)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(VariantAnnotation)
gwasvcf::set_bcftools("/share/home/zhangjunyu/Software/bcftools/bin/bcftools")
gwasvcf::set_plink()

outcome_list = fread("/share/home/linmiao/bio/rawdata/summary_stats_R10_manifest.tsv")
Files = list.files("/share/home/linmiao/bio/rawdata/UKB_PP_SNP/UKB_PP_SNP_1e5", pattern = ".vcf")
Files = gsub(".vcf", "", Files)
Files = Files

# 获取系统中的核心数并设置要使用的核心数量
num_cores <- detectCores()
num_cores_to_use <- min(30, num_cores)

# 定义并行任务函数
run_task <- function(id) {
  tryCatch({
    cat("Processing ID: ", id, "\n")
    expd <- gwasvcf::query_gwas(paste("/share/home/linmiao/bio/rawdata/UKB_PP_SNP/UKB_PP_SNP_1e5/", id, ".vcf", sep = ""), pval = 5e-8)
    expd <- gwasglue::gwasvcf_to_TwoSampleMR(expd, type = "exposure")
    expd <- expd %>% filter(eaf.exposure > 0.01)
    retain_snps <- expd %>%
      dplyr::select(rsid = SNP, pval = pval.exposure) %>%
      ieugwasr::ld_clump(., plink_bin = genetics.binaRies::get_plink_binary(), bfile = "/share/home/zhangjunyu/Rawdata/bfile/data_maf0.01_rs_ref") %>%
      {.$rsid}
    expd <- subset(expd, SNP %in% retain_snps) 
    
    cat("Finished processing exposure data for ID: ", id, "\n")

    outd <- gwasvcf::query_gwas(paste("/share/home/linmiao/bio/rawdata/Finngen/phenotype_vcf/finngen_R10_", outcome_id, ".vcf", sep = ""), rsid = expd$SNP, proxies = "yes", 
                                bfile = "/share/home/zhangjunyu/Rawdata/bfile/data_maf0.01_rs_ref")
    outd <- gwasglue::gwasvcf_to_TwoSampleMR(outd, "outcome")
    dat <- TwoSampleMR::harmonise_data(expd, outd)
    MR_result <- TwoSampleMR::mr(dat)

    cat("Finished MR analysis for ID: ", id, "\n")

    MR_result$id.exposure <- paste("UKB_Olink", gsub("_.*", "", id), sep = "_")
    MR_result$id.outcome <- outcome_id
    MR_result$exposure <- gsub("_.*", "", id)
    MR_result$outcome <- outcome
    MR_result <- MR_result[, c(1:2, 4, 3, 5:9)]

    write.csv(MR_result, paste(output_dir, "/", gsub("_.*", "", id), "_vs_", outcome_id, ".csv", sep = ""), row.names = FALSE)

    cat("Finished writing CSV for ID: ", id, "\n")
  }, error = function(e) {
    cat("Error encountered for ID ", id, ": ", e$message, "\n")
  })
}

# 遍历所有的 outcome_id 并行计算
for (outcome_id in unlist(outcome_list[, 1])[1:100]) {

  outcome = outcome_list[which(unlist(outcome_list[, 1]) == outcome_id), 2]
  outcome = gsub(",| |-|, |:|/|: |\\.|\\. ", "_", outcome)
  output_dir <- paste("/share/home/linmiao/bio/2024_UKB_pp/Drugtarget_D/UKB_Olink_FinnGen_", outcome_id, sep = "")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  }

  # 创建集群并指定日志输出文件
  cl <- makeCluster(num_cores_to_use, outfile = paste("/share/home/linmiao/bio/2024_UKB_pp/Drugtarget_D/Submitted_task/output/logfile", outcome_id, ".txt", sep = ""))

  # 在集群上运行初始化代码
  clusterEvalQ(cl, {
    library(genetics.binaRies)
    library(ieugwasr)
    library(gwasglue)
    library(gwasvcf)
    library(tidyverse)
    library(data.table)
    library(TwoSampleMR)
    library(VariantAnnotation)
    gwasvcf::set_bcftools("/share/home/zhangjunyu/Software/bcftools/bin/bcftools")
    gwasvcf::set_plink()
  })

  # 导出必要的变量和函数到集群节点
  clusterExport(cl, c("run_task", "outcome_list", "outcome_id", "outcome", "output_dir"))
  
  # 并行运行任务
  parLapply(cl, Files, run_task)
  
  # 停止集群
  stopCluster(cl)
}


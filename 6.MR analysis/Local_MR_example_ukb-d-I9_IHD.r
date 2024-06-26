library(genetics.binaRies)
library(ieugwasr)
library(gwasglue)
library(gwasvcf)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(VariantAnnotation)

Outcome_id = "ukb-d-I9_IHD" # 在这里改变结局选用的vcf文件
Outcome = "UKB_Ischaemic_heart_disease_wide_definition"

gwasvcf::set_bcftools("/share/home/zhangjunyu/Software/bcftools/bin/bcftools")
gwasvcf::set_plink()

Files = list.files("/share/home/linmiao/bio/rawdata/Finngen/Olink_vcf",pattern=".vcf")
Files = gsub(".vcf", "", Files)

for(id in Files){
    tryCatch({expd <- gwasvcf::query_gwas(paste("/share/home/linmiao/bio/rawdata/Finngen/Olink_vcf/", id, ".vcf", sep = ""), pval=5e-8)
    expd <- gwasglue::gwasvcf_to_TwoSampleMR(expd, type="exposure")
    retain_snps <- expd %>% dplyr::select(rsid=SNP, pval=pval.exposure) %>%
          ieugwasr::ld_clump(., plink_bin=genetics.binaRies::get_plink_binary(), bfile="/share/home/zhangjunyu/Rawdata/bfile/data_maf0.01_rs_ref") %>%
          {.$rsid}
    expd <- subset(expd, SNP %in% retain_snps)
    outd <- gwasvcf::query_gwas(paste("/share/home/linmiao/bio/rawdata/MR/IEU_GWAS/", Outcome_id, ".vcf.gz", sep = ""), rsid = expd$SNP, proxies="yes", 
                                                                                   bfile="/share/home/zhangjunyu/Rawdata/bfile/data_maf0.01_rs_ref")
    outd <- gwasglue::gwasvcf_to_TwoSampleMR(outd, "outcome")
    dat <- TwoSampleMR::harmonise_data(expd, outd)
    MR_result <- TwoSampleMR::mr(dat)
    MR_result$id.exposure <- paste("Finngen_Olink", gsub(".*_","",id), sep = "_")
    MR_result$id.outcome <- Outcome_id
    MR_result$exposure <- gsub(".*_","",id)
    MR_result$outcome <- Outcome
    MR_result <- MR_result[,c(1:2,4,3,5:9)]
    write.csv(MR_result, paste("/share/home/zhangjunyu/Project/IHD_Proteomic_analysis/Result/MR/Finngen_Olink_",Outcome,"/", gsub(".*_","",id), "_vs_", Outcome, ".csv",sep = ""), row.names = F)
},error=function(e){message("Error encountered for ID ", id)})
}

#预处理
#需要把vcf文件变成vcf.gz文件，并生成index文件
#bgzip finngen_R10_C3_BLADDER_EXALLC.vcf
#需要对序列排序，bcftools sort -Oz -o sorted_finngen_R10_C3_BLADDER_EXALLC.vcf.gz finngen_R10_C3_BLADDER_EXALLC.vcf.gz
#然后创建索引 bcftools index sorted_finngen_R10_C3_BLADDER_EXALLC.vcf.gz

library(gwasglue) 
library(tidyverse) 
library(coloc)
library(geni.plots)
library(data.table)

gwasvcf::set_bcftools("/share/home/zhangjunyu/Software/bcftools/bin/bcftools")

coloc_to_gassocplot_geni.plots <- function(coloclist)
{
  markers <- dplyr::tibble(
    marker = coloclist$dataset1$snp,
    chr = coloclist$dataset1$chr,
    pos = coloclist$dataset1$pos,
    z_1 = coloclist$dataset1$z,
    z_2 = coloclist$dataset2$z
  )
  message("Extracting LD matrix for ", nrow(markers), " variants")
  ld <- ieugwasr::ld_matrix(markers[["marker"]], with_alleles=FALSE, 
                bfile="/share/home/zhangjunyu/Rawdata/bfile2/EUR", plink_bin = genetics.binaRies::get_plink_binary())
  message("Found ", nrow(ld), " variants in LD reference panel")
  index <- match(rownames(ld), markers[["marker"]])
  markers <- markers[index, ]
  stopifnot(all(markers$marker == rownames(ld)))
  list(markers = markers, corr = ld) %>% return()
}

Outcome_id = "finngen_R9_DM_NEPHROPATHY_EXMORE" # 在这里改变结局选用的vcf文件
Outcome = "Diabetic_nephropathy"
outcomne_info <- fread("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/conc_vcf/R9_manifest.tsv")
outcome_info <- outcomne_info %>% filter(phenocode == gsub("finngen_R9_", "", Outcome_id))

Files = list.files("/share/home/linmiao/bio/rawdata/UKB_PP_SNP/DataWashvcf/Data",pattern = "^ADM_|^AGER_|^ANGPTL3_|^CCER2_|^CD160_|^CD300A_|^CFC1_|^CFD_|^CGA_|^CHRDL1_|^CX3CL1_|^DCBLD2_|^FAM3C_|^FOLR1_|^FSTL1_|^GABARAP_|^GCNT1_|^GPR37_|^HIP1R_|^HLA-E_|^IFNGR1_|^IGFBP1_|^IGFBP4_|^MYOC_|^NBL1_|^NCR1_|^NTproBNP_|^PALM2_|^POLR2F_|^PTGDS_|^RNASET2_|^SCARF1_|^SCGB3A2_|^SHISA5_|^TNFRSF1B_|^TNFRSF4_|^VWC2L_")
Files = gsub(".vcf", "", Files)

protein_list <- fread("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/conc_vcf/gencode.v37lift37.annotation.gtf")
protein_list <- protein_list %>% filter(V3 == "gene")

coloc_result <- tibble()

for (id in Files){
  tryCatch({
  protein_id <- gsub("_.*", "", id)
  #protein_GWAS <- fread(paste("/share/home/linmiao/bio/rawdata/UKB_PP_SNP/DataWashvcf/Data/",id, ".vcf", sep=""))
  #protein_GWAS$Study_1 <- gsub(".*:.*:.*:(.*):.*:.*", "\\1", protein_GWAS$Study_1)
  #protein_GWAS$Study_1 <- as.numeric(protein_GWAS$Study_1)
  #max_row <- protein_GWAS[which.max(Study_1)]
  #names(max_row)[1] <- "CHROM"
  #chrpos <- paste(max_row$CHROM, ":", max_row$POS-90000, "-", max_row$POS+90000, sep="")
  if (protein_id == "NTproBNP"){
    pro_loc <- protein_list[grep(paste("gene_name \"", "NPPB", "\"" ,sep=""),V9),]
  }else if (protein_id == "PALM2") {
    pro_loc <- protein_list[grep(paste("gene_name \"", "PALM2AKAP2", "\"" ,sep=""),V9),]
  }else {pro_loc <- protein_list[grep(paste("gene_name \"", protein_id, "\"" ,sep=""),V9),]}

  gen_mid <- floor((pro_loc$V4[1]+pro_loc$V5[1])/2)
  chrpos <- paste(gsub("chr","",pro_loc[1,1]), ":", gen_mid-90000, "-", gen_mid+90000, sep="")
  chrpos <- gsub("X", "23", chrpos)

  GWAS1 <- paste("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/conc_vcf/proteins_vcf/sorted_", id, ".vcf.gz",  sep = "")
  type1 <- "quant" #定量性状采用"quant"，定性形状采用“cc”
  GWAS2 <- paste("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Data/conc_vcf/DKD/sorted_", Outcome_id, ".vcf.gz",  sep = "")
  type2 <- "cc"    #定量性状采用"quant"，定性形状采用“cc”

  o <- gwasvcf::vcflist_overlaps(list(GWAS1, GWAS2), chrpos)

  GWAS1 <- o[[1]]
  GWAS2 <- o[[2]]
  if(length(GWAS1) == 0 | length(GWAS2) == 0) stop("No overlaps in GWAS document")
  if(length(GWAS1) != length(GWAS2)) stop("Refined GWAS has different length")
  tab1 <- GWAS1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
  tab2 <- GWAS2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
  index <- as.character(tab1$REF) == as.character(tab2$REF) &
        as.character(tab1$ALT) == as.character(tab2$ALT) &
        as.character(tab1$seqnames) == as.character(tab2$seqnames) &
        tab1$start == tab2$start & tab1$ID!="<NA>" & tab2$ID!="<NA>"
  stopifnot(sum(index) > 0)
  tab1$AF[is.na(tab1$AF)] <- 0.5
  tab2$AF[is.na(tab2$AF)] <- 0.5
  #tab1$SS <- #如果sample size有缺失，需要补全
  tab2$SS <- outcome_info$num_controls + outcome_info$num_cases  #如果sample size有缺失，需要补全
  tab2$NC <- outcome_info$num_cases  #如果NC(number of cases)有缺失，需要补全
  out1 <- tab1[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, 
                                type = type1, snp = names(GWAS1)[index], z = .$ES / .$SE, chr = .$seqnames, 
                                pos = .$start)}
  out2 <- tab2[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, 
                                type = type2, snp = names(GWAS2)[index], z = .$ES / .$SE, chr = .$seqnames, 
                                pos = .$start)}
  if(type1 == "cc")	out1$s <- mean(tab1$NC / tab1$SS, na.rm=TRUE)
  if(type2 == "cc") out2$s <- mean(tab2$NC / tab2$SS, na.rm=TRUE)

  #进行共定位分析               
  vres <- coloc::coloc.abf(out1, out2)
  coloc_result <- bind_rows(coloc_result, vres$summary)

  #可视化1：两个GWAS的region合并作图，会把两个trait的放在一张图上面
  temp <- coloc_to_gassocplot_geni.plots(list(dataset1=out1, dataset2=out2))
  p <- fig_region_stack(
    data = temp$markers,
    traits = c(protein_id, gsub("_", "", Outcome)), #在这里可以更改可视化图片的title
    corr = temp$corr,
    build = 37,
    # highlights = "rs11265611",        #可以高亮想要高亮的的SNP
    title_center = TRUE
  )
  ggsave(paste("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/coloc/", Outcome_id, "/",protein_id,"_", Outcome,".pdf",sep=""), plot = p, device = "pdf")

},error= function(e) {
    print(paste("Caught an error:", e$message))
  })
}
#coloc_result$Proteins = gsub("_.*", "", Files)[-10]
#coloc_result$Outcome = Outcome
#coloc_result <- coloc_result[, c(7:8, 1:6)]
#write.csv(coloc_result, paste("/share/home/zhangjunyu/Project/Diabetes_Proteomic_analysis/Result/conc/", Outcome_id, "/",Outcome,".csv",sep=""), row.names = F)


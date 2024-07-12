#预处理
#需要把vcf文件变成vcf.gz文件，并生成index文件
#bgzip finngen_R10_C3_BLADDER_EXALLC.vcf
#需要对序列排序，bcftools sort -Oz -o sorted_finngen_R10_C3_BLADDER_EXALLC.vcf.gz finngen_R10_C3_BLADDER_EXALLC.vcf.gz
#然后创建索引 bcftools index sorted_finngen_R10_C3_BLADDER_EXALLC.vcf.gz

library(gwasglue) 
library(dplyr) 
library(coloc)
library(geni.plots)

gwasvcf::set_bcftools("/share/home/zhangjunyu/Software/bcftools/bin/bcftools")

#设置比对的范围
chrpos <- "19:58751560-58931560"

#设置进行比对的两个GWAS文件
GWAS1 <- "/share/home/zhangjunyu/Practice/240711_conc/A1BG_P04217_OID30771_v1_Inflammation_II.vcf.gz"
type1 <- "quant" #定量性状采用"quant"，定性形状采用“cc”
GWAS2 <- "/share/home/zhangjunyu/Practice/240711_conc/sorted_finngen_R10_C3_BLADDER_EXALLC.vcf.gz"
type2 <- "cc"    #定量性状采用"quant"，定性形状采用“cc”

#读取两个GWAS文件，并将其洗为conc要求格式的文件
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
			tab1$start == tab2$start
stopifnot(sum(index) > 0)
tab1$AF[is.na(tab1$AF)] <- 0.5
tab2$AF[is.na(tab2$AF)] <- 0.5
#tab1$SS <- #如果sample size有缺失，需要补全
tab2$SS <- 412181  #如果sample size有缺失，需要补全
tab2$NC <- 2193  #如果NC(number of cases)有缺失，需要补全
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
#如果确实找到了共定位H4的有力证据，我们可以提取每个SNP以H4为因果条件的后验概率.可以通过以下方式提取更可能的因果变体
subset(vres$results,SNP.PP.H4>0.01)
#95% 可信度设定
o <- order(vres$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(vres$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
vres$results[o,][1:w,]$snp

#可视化1：两个GWAS的region合并作图，会把两个trait的放在一张图上面
#参考：https://github.com/MRCIEU/gwasglue/blob/c2d5660eed389e1a9b3e04406b88731d642243f1/R/gassocplot.r
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

temp <- coloc_to_gassocplot_geni.plots(list(dataset1=out1, dataset2=out2))
pdf("conc_show1.pdf")
fig_region_stack(
  data = temp$markers,
  traits = c("A1BG", "BLADDER_EXALLC"), #在这里可以更改可视化图片的title
  corr = temp$corr,
  build = 37,
  # highlights = "rs11265611",        #可以高亮想要高亮的的SNP
  title_center = TRUE
)
dev.off()

#可视化2：单一region作图,可能大部分时候并不能用到，如果想得到两个子图，自行拼接，可以考虑这个
#参考：https://github.com/MRCIEU/gwasglue/blob/c2d5660eed389e1a9b3e04406b88731d642243f1/R/gassocplot.r
coloc_to_gassocplot_geni.plots <- function(coloclist)
{
  markers <- dplyr::tibble(
    marker = coloclist$dataset1$snp,
    chr = coloclist$dataset1$chr,
    pos = coloclist$dataset1$pos
  )
  z <- dplyr::tibble(
    id1 = coloclist$dataset1$z,
    id2 = coloclist$dataset2$z
  )
  message("Extracting LD matrix for ", nrow(markers), " variants")
  ld <- ieugwasr::ld_matrix(markers[["marker"]], with_alleles=FALSE, 
                bfile="/share/home/zhangjunyu/Rawdata/bfile2/EUR", plink_bin = genetics.binaRies::get_plink_binary())
  message("Found ", nrow(ld), " variants in LD reference panel")
  index <- match(rownames(ld), markers[["marker"]])
  markers <- markers[index, ]
  z <- z[index, ]
  stopifnot(all(markers$marker == rownames(ld)))
  list1 = list(markers = bind_cols(markers,z[,1]), corr = ld)
  list1$markers <- rename(list1$markers, 'z'= 'id1')
  list2 = list(markers = bind_cols(markers,z[,2]), corr = ld)
  list2$markers <- rename(list2$markers, 'z'= 'id2')
  result = list(list1, list2)
  return(result)
}

temp <- coloc_to_gassocplot_geni.plots(list(dataset1=out1, dataset2=out2))
pdf("conc_show2.pdf")
fig_region(
  data = temp[[1]]$markers,
  corr = temp[[1]]$corr,
  build = 37
)
fig_region(
  data = temp[[2]]$markers,
  corr = temp[[2]]$corr,
  build = 37
)
dev.off()



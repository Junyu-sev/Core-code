import gwaslab as gl
import pandas as pd

lst = pd.read_csv("/share/home/linmiao/bio/rawdata/UKB_PP_SNP/DataWashtsv/1Oncology-Filelist.tsv")
lst = lst.to_numpy()

for gwas_id in lst:
    mysumstats = gl.Sumstats('/share/home/linmiao/bio/rawdata/UKB_PP_SNP/DataWashtsv/1Oncology/' + gwas_id[0] + '.tsv', 
             snpid="snpid",
             chrom="CHROM",
             pos="GENPOS",
             ea="ALLELE1",
             nea="ALLELE0",
             eaf="A1FREQ",
             beta="BETA",
             se="SE",
             mlog10p="LOG10P",
             n="N",
             build="38")
    mysumstats.basic_check()
    mysumstats.assign_rsid( n_cores = 35,
                        ref_rsid_vcf ="/share/home/zhangjunyu/Rawdata/Reference_dbSNP/GCF_000001405.40.gz",
                        chr_dict = gl.get_number_to_NC(build="38"))
    mysumstats.liftover(n_cores=35, 
                    from_build="38", 
                    to_build="19",
                    remove=True)
    mysumstats.to_format('/share/home/linmiao/bio/rawdata/UKB_PP_SNP/DataWashvcf/1Oncology/'+gwas_id[0], fmt="vcf", build="19",id_use="rsID",xymt_number=True)#可以再研究一下，以什么格式输出，后面在研究一下



# Data Clean for FinnGen GWAS Summary

## Background

**[TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/)** package is used to perform **two-sample Mendelian Randomization**. Before local analysis, GWAS summary data needs to be adjusted to standard form. **[IEU website](https://gwas.mrcieu.ac.uk/)** offers various standarlized GWAS summary as the form of VCF documents. However, FinnGen GWAS summary data is not provided and here is a tutorial for FinnGen data clean.

## About FinnGen

**[FinnGen](https://www.finngen.fi/en/access_results)** is a large public-private partnership aiming to collect and analyse genome and health data from 500,000 Finnish biobank participants. FinnGen results can be browsed and summary statistics data can be downloaded from Google cloud storage free of charge.The latest data is released in December 2023. Latest data profile is as follows:

- Total sample size: **412,181** (230,310 females and 181,871 males)
- Total number of variants analyzed: **21,311,942 variants**
- Number of disease endpoints (phenotypes) available: **​2,408 endpoints**
- LOF burden test results,coding variant results and **proteomics QTL** results are also included

More detailed description for FinnGen data could be found [here](https://finngen.gitbook.io/documentation). **All accessible FinnGen data were downloaded and could be found in `/share/home/linmiao/bio/rawdata/Finngen`.** And index for FinnGen data could be download in below three linkages:

1. [GWAS summary index for various endpoints](https://storage.googleapis.com/finngen-public-data-r10/summary_stats/R10_manifest.tsv)
2. [pQTL index from Olink proteomics](https://console.cloud.google.com/storage/browser/_details/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/probe_map.tsv?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)))
3. [pQTL index from Somascan proteomics](https://console.cloud.google.com/storage/browser/_details/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/probe_map.tsv?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)))

## Tools for data clean

**[GWASlab](https://cloufield.github.io/gwaslab/)** is a handy Python toolkit for handling GWAS summary statistics (sumstats). We will use it to **fix lost SNPID information, liftover positions of all variants from hg38 to hg19, and generate VCF documents** for every FinnGen GWAS summary.

To begin with, you need to install this toolkit. Please follow the tutorials from **[GWASlab website](https://cloufield.github.io/gwaslab/#install)**.

Then, please enter the path where you install GWASlab and find file `data/formatbook.json`. Just for example, my path is `/share/home/zhangjunyu/anaconda3/envs/gwaslab/lib/python3.9/site-packages/gwaslab`). You need to add one line of code to this file as below (this is to guarantee the correct form of VCF to be output).

```java
    ......
    "vcf": {
        "meta_data": {
            "format_name": "vcf",
            "format_source": "https://github.com/MRCIEU/gwas-vcf-specification/tree/1.0.0",
            "format_version": 20220923,
            "format_citation": "Lyon, M.S., Andrews, S.J., Elsworth, B. et al. The variant call format provides efficient and robust storage of GWAS summary statistics. Genome Biol 22, 32 (2021). https://doi.org/10.1186/s13059-020-02248-0",
            "format_fixed_header": "##fileformat=VCFv4.2\r\n##FILTER=<ID=PASS,Description=\"All filters passed\">\r\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\r\n##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size estimate relative to the alternative allele\">\r\n##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error of effect size estimate\">\r\n##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 p-value for effect estimate\">\r\n##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency in the association study\">\r\n##FORMAT=<ID=SS,Number=A,Type=Float,Description=\"Sample size used to estimate genetic effect\">\r\n##FORMAT=<ID=EZ,Number=A,Type=Float,Description=\"Z-score provided if it was used to derive the EFFECT and SE fields\">\r\n##FORMAT=<ID=SI,Number=A,Type=Float,Description=\"Accuracy score of summary data imputation\">\r\n##FORMAT=<ID=NC,Number=A,Type=Float,Description=\"Number of cases used to estimate genetic effect\">\r\n##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Study variant identifier\">\r\n##META=<ID=TotalVariants,Number=1,Type=Integer,Description=\"Total number of variants in input\">\r\n##META=<ID=VariantsNotRead,Number=1,Type=Integer,Description=\"Number of variants that could not be read\">\r\n##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description=\"Total number of harmonised variants\">\r\n##META=<ID=VariantsNotHarmonised,Number=1,Type=Integer,Description=\"Total number of variants that could not be harmonised\">\r\n##META=<ID=SwitchedAlleles,Number=1,Type=Integer,Description=\"Total number of variants strand switched\">\r\n##META=<ID=TotalControls,Number=1,Type=Integer,Description=\"Total number of controls in the association study\">\r\n##META=<ID=TotalCases,Number=1,Type=Integer,Description=\"Total number of cases in the association study\">\r\n##META=<ID=StudyType,Number=1,Type=String,Description=\"Type of GWAS study [Continuous or CaseControl]\">",
            "format_fixed": [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT"
            ],
            "format_format": [
                "ID",
                "AF",     #### Left Is The Only Thing to Add, Dont Add This Annotation ! ####
                "SS",
                "ES",
                "SE",
                "LP",
                "SI",
                "EZ"
            ],
    ......
```

Also, you need to download **reference dbSNP** which will be utilized to harmaonize GWAS summary data and fix lost rsID. The most recently published dbSNP file could be found in [NCBI](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/). You can download it by:

```linux
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
```
This is quit a time-spending matter. You can also download this file in your local environment and upload it into HPC rather than directly downloading it in HPC (maybe less time-costing).

Then you need to make index for this reference file. If you have installed **bcftools** in your environment, you can use the code `bcftools index /Your/Path/To/Reference_dbSNP/GCF_000001405.40.gz`. If you haven't install bcftools, you can finish installation following [tutorial here](http://www.htslib.org/download/).

All preparation is made, time to clean FinnGen data. **Activate your conda environment and python interpretor.** Following example codes are all you need to make data clean. You could get [code instructions here](https://cloufield.github.io/gwaslab/AssignrsID/). Additional changes could be made according to your actual situation.

```python
import gwaslab as gl

mysumstats = gl.Sumstats("/Your/FinnGen/GWAS/Path.gz", 
        snpid="rsids",
        chrom="#chrom",
        pos="pos",
        ea="alt",
        nea="ref",
        eaf="af_alt",
        beta="beta",
        se="sebeta",
        p="pval",
        mlog10p="mlogp",
        build="38")
mysumstats.basic_check()
mysumstats.assign_rsid(n_cores = 20,
                        ref_rsid_vcf ="/Your/Path/to/Reference_dbSNP/GCF_000001405.40.gz",
                        chr_dict = gl.get_number_to_NC(build="38"))
mysumstats.liftover(n_cores=20, 
                    from_build="38", 
                    to_build="19",
                    remove=True)
mysumstats.to_format("/Where/To/Save/Cleaned/Data/filename", fmt="vcf", build="19",id_use="rsID",xymt_number=True)
```

## Finally

Some mistakes may exist above. Feel it free to feed back!

## 补充更新

gwaslab读取.vcf文件会报错，需要修改/share/home/zhangjunyu/anaconda3/envs/gwaslab/lib/python3.9/site-packages/gwaslab/io_preformat_input.py文件中get_readargs_header函数的定义，改成下面这个就可以啦~ 读取.vcf.gz文件不受任何影响

```python
def get_readargs_header(inpath,readargs):
    if "vcf.gz" in inpath:
        with gzip.open(inpath,'r') as file:      
            skip=0
            for line in file:        
                if line.decode('utf-8').startswith('##'):
                    skip+=1
                else:
                    readargs["skiprows"]=skip
                    readargs["sep"]="\t"
                    break
    elif "vcf" in inpath:
        with open(inpath, 'r') as file:
            skip = 0
            for line in file:
                if line.startswith('##'):
                    skip += 1
                else:
                    readargs["skiprows"] = skip
                    readargs["sep"] = "\t"
                    break
    readargs_header = readargs.copy()
    readargs_header["nrows"]=1
    readargs_header["dtype"]="string"
    return readargs_header
```

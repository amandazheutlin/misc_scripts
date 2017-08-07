# gene set analysis of change phenotype
# calcium genes (N = 26)

#### housekeeping
set.seed(1)
workdir <- "/data/swe_gwas/ABZ/NAPLS/CAgenes/"
genedir <- "/data/swe_gwas/ABZ/NAPLS/postqc/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2", "SKAT")
invisible(lapply(libs, require, character.only = TRUE))

#### get calcium genes from UCSC
# group = Genes and Gene Predictions
# track = GENCODE v22; table = kgXref
# paste list of identifiers (gene symbols)
# output format = selected fields
# get output -> choose linked table 'knownCanonical'
# check chrom, chromStart, chromEnd -> get output

# CA gene list 
# CA         <- read.table("CAgenes1.txt",header=F)
# CApos        <- read.table("CAgenes_pos.txt",header=T)
# names(CApos) <- c("UCSC","symbol","chr","start","end")
# CApos        <- CApos[!(CApos$chr=="n/a"),]
# CApos$chrno  <- substr(CApos$chr,4,length(CApos$chr))
# write.table(CApos[,c(6,4,5,2)],"CApos_plink.txt",sep="\t",
#             col.names=F,row.names=F,quote=F)

##################################### SKAT tests
# let's try SKAT to test for common and rare variant burden in CA genes
# steps: 1) generate SNP file (SSD), 2) get genotype data from SSD
# 3) read in phenotype and covariates, 4) SKAT for common + rare

#### generate SSD + open it
# set ID
# map          <- read.table("CAsnps-exome.map",header=F)
# setID        <- map[,2,drop=F]
# setID$SetID  <- "CA"
# setID        <- setID[,c(2,1)] # set ID, SNP ID
# 
# write.table(setID,"CA_setID.txt",col.names=F,row.names=F,quote=F,sep="\t")

Generate_SSD_SetID(paste0(genedir,"NAPLS_zcall_011013_10.bed"),
                   paste0(genedir,"NAPLS_zcall_011013_10.bim"),
                   paste0(genedir,"NAPLS_zcall_011013_10.fam"),
                   paste0(workdir,"CA_setID.txt"),
                   "SSDfile.CAgenes",
                   "SSDfile.CAgenes.info")

SSD <- Open_SSD("SSDfile.CAgenes","SSDfile.CAgenes.info")
# Close_SSD()

#### get genotypes from SSD
# 1 is for set 1 (only 1 set)
geno <- Get_Genotypes_SSD(SSD, 1, is_ID=T)

#### read in phenotype and covar files

# add phenotype to fam file
# fam           <- read.table(paste0(genedir,"NAPLS_zcall_011013_10.fam"),header=F)
# names(fam)    <- c("FID","IID","PID","MID","Sex","Phenotype")
# gm            <- read.table("../change-pheno.txt",header=F)
# names(gm)     <- c("FID","IID","Phenotype")
# fam$Phenotype <- gm[match(fam$FID,gm$FID),3]

# write.table(fam,paste0(genedir,"NAPLS_zcall_011013_10_pheno.fam"),col.names=F,
#             row.names=F,quote=F,sep="\t")

# this reads in FID (actual ID), IID ('1' for everyone), PID, MID, sex,
# GM pheno, age, sex (lol), 5 PCs
phenocov <- Read_Plink_FAM_Cov(paste0(genedir,"NAPLS_zcall_011013_10_pheno.fam"),
                   "../covar-file.txt", Is.binary=F, cov_header=F)
pheno    <- phenocov[,6]
cov      <- phenocov[,7:13] %>% as.matrix()

#### SNP-set kernal association test (SKAT) test of common/rare burden

# the null model
# formula = y.b ~ X
nullmodel <- SKAT_Null_Model(pheno ~ cov,out_type="C",n.Resampling = 1000,
                             type.Resampling = "bootstrap", Adjustment = F)

# association test
# geno = matrix of SNPs; nullmodel = the null model
# method is the type of sum test (C = combined)
# rare / common beta weights are set to default
# r.corr parameters set to 0 (default) for SKAT (1 is burden)
# rare vs. common cutoff -> above is common, below is rare
assoc   <- SKAT_CommonRare(geno, nullmodel, method="C", test.type="Joint",
                           CommonRare_Cutoff = .05)
burden  <- SKAT_CommonRare(geno, nullmodel, method="C", test.type="Joint",
                           CommonRare_Cutoff = .05, 
                           r.corr.rare = 1, r.corr.common = 1)
common  <- SKAT_CommonRare(geno, nullmodel, method="C",test.type="Common.Only")
rare    <- SKAT_CommonRare(geno, nullmodel, method="C",test.type="Rare.Only")

assoc$p.value
burden$p.value
common$p.value
rare$p.value




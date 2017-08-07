# exome x imaging change phenotype
# NAPLS subjects with imaging + exome (N = 345)

# housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2")
invisible(lapply(libs, require, character.only = TRUE))

# covar file
covar <- read.table(paste0(workdir,"age-sex.txt"),header=T)
mds <- read.table(paste0(workdir,"NAPLS-mds.mds"),header=T,stringsAsFactors = F)

covar2 <- merge(covar,mds[,c(1,4:8)],by="FID")
write.table(covar2,"covar-file.txt",col.names=F,row.names=F,sep="\t",quote=F)

#####################################################
# association between exome and change pheno
# 238,689 variants; 345 subjects
# covar = age, sex, 5 PCs
assoc <- read.table("exome_changepheno_allsubs.assoc.linear",header=T)
assoc <- assoc[order(assoc$P),]
assoc05 <- subset(assoc,P<.05)
assoc05$P.bon <- p.adjust(assoc05$P,method="bon",n=238689)
assocbon <- subset(assoc05,P.bon<.05) # 607 markers

snps <- assocbon[,2,drop=F]
write.table(snps,"bon-snps_changeassoc.txt",col.names=F,
            row.names=F,quote=F,sep="\t")

# map exome markers to chr positions
map <- read.table(paste0(workdir,"postqc/NAPLS_zcall_011013_10.bim"),header=F)
names(map) <- c("chr","snp","dis","pos","A1","A2")

positions <- assocbon[,1:2]
positions$CHR2 <- paste0("chr",positions$CHR)
positions$POS <- map[match(assocbon$SNP,map$snp),4]
positions$POS2 <- positions$POS + 1

write.table(positions[,c(3:5)],"snp-bon_positions.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

# errors from UCSC (start pos)
errors <- c(145741938,146156729,139276351,135083881,
            58772640,146278811,47836593,50943916,15905,
            139751649,145140216,51049171,9966)
pos2 <- positions[!(positions$POS %in% errors),]

write.table(pos2[,c(3:5)],"snp-bon_positions_noerrors.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

# read in UCSC output
# markers annotated with UCSC gene names
genes <- read.table("ucsc_output.txt",header=F)
genes <- genes[,c(1:5,7,10,13)] # clean output
names(genes) <- c("chr","GENCODE","region","pos1","pos2",
                  "strand.probably","UCSCgene","UCSCtranscript")

genelist <- unique(genes$UCSCgene)
write.table(genelist,"UCSCgenes.csv",sep=",",quote=F,col.names=F,row.names=F)

# LD markers
LD <- read.table("bonsnps-LD.prune.in",header=F)
names(LD) <- "snp"
LD$CHR <- paste0("chr",map[match(LD$snp,map$snp),1])
LD$POS <- map[match(LD$snp,map$snp),4]
LD$POS2 <- LD$POS + 1

LD2 <- LD[!(LD$POS %in% errors),] # don't get duplicate errors, at least

write.table(LD2[,c(2:4)],"snp-bon_positions_LD.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

# ucsc output - LD
genesld <- read.table("ucsc_output_ld.txt",header=F)
names(genesld)[1:4] <- c("chr","pos1","pos2","gene")

genelist2 <- unique(genesld$gene)
write.table(genelist2,"UCSCgenes_ld.csv",sep=",",quote=F,col.names=F,row.names=F)

# ucsc output - LD - gene symbols
sym <- read.table("ucsc_output_ld_sym.txt",header=F)
names(sym) <- c("ucsc_gene","gene_symbol")
genesld$gene_symbol <- sym[match(genesld$gene,sym$ucsc_gene),2] %>% as.character

genelist3 <- unique(genesld$gene_symbol) 
write.table(genelist3,"genesymbols_ld.csv",sep=",",quote=F,col.names=F,row.names=F)

genes$gene_symbol <- sym[match(genes$UCSCgene,sym$ucsc_gene),2] %>% as.character
genelist4 <- unique(genes$gene_symbol) 
write.table(genelist4,"genesymbols.csv",sep=",",quote=F,col.names=F,row.names=F)

# allele frequencies of 607 bon snps
freq <- read.table("bonsnps-cp-A1freq.frq",header=T)
freq <- freq[order(freq$MAF),]

ggplot(freq,aes(x=MAF)) + 
  geom_histogram(stat="bin",binwidth = .001, position="dodge")

#####################################################
# association between exome and change pheno
# 42,035 variants (MAF > .01); 345 subjects
# covar = age, sex, 5 PCs
assoc <- read.table("exome_gwas/exome_cp_common.assoc.linear",header=T)
assoc <- assoc[order(assoc$P),]
assoc05 <- subset(assoc,P<.05)
assoc05$P.bon <- p.adjust(assoc05$P,method="bon",n=42035)
assocbon <- subset(assoc05,P.bon<.05) # 28 markers

# map exome markers to chr positions
map <- read.table(paste0(workdir,"postqc/NAPLS_zcall_011013_10.bim"),header=F)
names(map) <- c("chr","snp","dis","pos","A1","A2")

snps <- assocbon[,2,drop=F]
write.table(snps,"bon-commonsnps_changeassoc.txt",col.names=F,
            row.names=F,quote=F,sep="\t")

pos <- assocbon[,1:2]
pos$CHR <- paste0("chr",pos$CHR)
pos$POS <- map[match(pos$SNP,map$snp),4]
pos$POS2 <- pos$POS + 1

write.table(pos[,c(1,3,4)],"commonsnp-bon_positions.txt",
            col.names=T,row.names=F,quote=F,sep="\t")

# converted to ucsc
common <- read.table("ucsc_output_common.txt",header=F)
names(common)[1:4] <- c("chr","pos1","pos2","gene")

# coverted to gene symbol
sym <- read.table("ucsc_output_ld_sym.txt",header=F)
names(sym) <- c("ucsc_gene","gene_symbol")
common$gene_symbol <- sym[match(common$gene,sym$ucsc_gene),2] %>% as.character

commonlist <- unique(common$gene_symbol) 
write.table(commonlist,"genesymbols_common.csv",sep=",",quote=F,col.names=F,row.names=F)

# covert ucsc to NM_ID
nm <- read.table("commonbonsnps_ucsc.txt",header=F)
write.table(nm[,4],"nm_commonbonsnps.txt",col.names=F,
            row.names=F,quote=F,sep="\t")

# baseline data for 'replication'
# rh_superiorfrontal_thickness
bl        <- read.csv("exome_gwas/NAPLS2_fs_bl_all_new_4_10_2016.csv",
                      header=T,sep="\t",stringsAsFactors = F)
covardisc <- read.table("covar-file.txt",header=F,stringsAsFactors = F)
discovery <- covardisc[,1]

rhsft     <- bl[,c("imagingID","rh_superiorfrontal_thickness")]
rhsft$IID <- rep(1,750)
names(rhsft)[1] <- "FID"
rhsft     <- rhsft[,c(1,3,2)]
repsample <- rhsft[!(rhsft$FID %in% discovery),] # N = 406

demo      <- read.csv("exome_gwas/NAPLS2_demo_bl_all_cleaned.csv",
                      header=T,sep="\t",stringsAsFactors = F)
repcovar  <- mds[!(mds$FID %in% discovery),] %>% dplyr::select(c(1:2,4:8))
repcovar  <- merge(repcovar,demo[,c(2,7,13,15)],by="FID")
repcovar$SiteInt <- as.factor(repcovar$SiteInt) %>% as.numeric() # covs must be numeric

write.table(repsample,"exome_gwas/bl_rep_pheno.txt",col.names=F,
            row.names=F,quote=F,sep="\t")
write.table(repcovar,"exome_gwas/bl_rep_covar.txt",col.names=F,
            row.names=F,quote=F,sep="\t")

# association
repassoc <- read.table("exome_gwas/combonsnps_bl_cov.assoc.linear",header=T)
repassoc <- repassoc[order(repassoc$P),]





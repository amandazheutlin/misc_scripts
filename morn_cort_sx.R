# SLEEP DYSREGULATION IN SCZ
# correlations between morningness + cortical thickness, sx
# also look at CLOCK genes + cortical thickness, sx
# swedish twins

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/circadian_rhythms/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2", "nlme")
invisible(lapply(libs, require, character.only = TRUE))

#### load data
morn <- read.table("morningness.txt",header=T)
cort <- read.table("/data/swe_gwas/ABZ/MRI_ROI_swedes/FS_VolMeanThick_20151120.txt",header=T)
cort$IID <- substr(cort$scanid,1,10)

morn$lthick     <- cort[match(morn$PTID,cort$IID),20]
morn$rthick     <- cort[match(morn$PTID,cort$IID),21]
morn$morn.score <- (6-morn$BestGetUp) + (6-morn$BestGoBed) +morn$HwEsyUp +morn$HwEsyBd + (5-morn$TypePeop)

#### correlations between morningness + cortical thickness
m1 <- lme(lthick ~ morn.score, random = ~1|FID, data = morn, na.action = na.omit)
m2 <- lme(rthick ~ morn.score, random = ~1|FID, data = morn, na.action = na.omit)

patients <- subset(morn,dx==1 | dx==3)
m3 <- lme(lthick ~ morn.score, random = ~1|FID, data = patients, na.action = na.omit)
m4 <- lme(rthick ~ morn.score, random = ~1|FID, data = patients, na.action = na.omit)

SZ <- subset(morn,dx==1)
m5 <- lme(lthick ~ morn.score, random = ~1|FID, data = SZ, na.action = na.omit)
m6 <- lme(rthick ~ morn.score, random = ~1|FID, data = SZ, na.action = na.omit)

#### correlations between thickness + morningness items
mornitems <- names(morn)[8:12]

# left thickness
lthick <- matrix(nrow=length(mornitems),ncol=4)
for (i in 1:length(mornitems)){
  f <- formula(paste("lthick", mornitems[i], sep="~"))
  m1 <- lme(fixed=f, random = ~1|FID, data = morn, na.action = na.omit)
  tval <- summary(m1)$tTable[2,4]
  lthick[i,1] <- mornitems[i]
  lthick[i,2] <- tval
  lthick[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  lthick[i,4] <- summary(m1)$tTable[2,5]
}

lthick        <- as.data.frame(lthick)
names(lthick) <- c("item","t","corr","p")
lthick$p.fdr  <- p.adjust(lthick$p, method="BH", n=5)

# right thickness
rthick <- matrix(nrow=length(mornitems),ncol=4)
for (i in 1:length(mornitems)){
  f <- formula(paste("rthick", mornitems[i],sep="~"))
  m1 <- lme(fixed=f, random = ~1|FID, data = morn, na.action = na.omit)
  tval <- summary(m1)$tTable[2,4]
  rthick[i,1] <- mornitems[i]
  rthick[i,2] <- tval
  rthick[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  rthick[i,4] <- summary(m1)$tTable[2,5]
}

rthick        <- as.data.frame(rthick)
names(rthick) <- c("item","t","corr","p")
rthick$p.fdr  <- p.adjust(rthick$p, method="BH", n=5)

#### correlations between clock SNPs and cortical thickness
# SNPs from GWAS: snp1 = rs11121022, snp2 = rs35833281
# 2 minor alleles = 2, hets = 1, majors = 0
clock <- read.table("/data/swe_gwas/ABZ/imputed_data/OMNI/call_genos/clocksnps.ped",header=F)
clock <- clock[,c(1:2,7:10)]
names(clock) <- c("FID","IID","snp1.A1","snp1.A2","snp2.A1","snp2.A2")

clock$snp1     <- paste0(clock$snp1.A1,clock$snp1.A2)
clock$snp2     <- paste0(clock$snp2.A1,clock$snp2.A2)
clock$snp1.var <- ifelse(clock$snp1=="11",2,ifelse(clock$snp1=="12",1,0))
clock$snp2.var <- ifelse(clock$snp2=="11",2,ifelse(clock$snp2=="12",1,0))

morn$snp1      <- clock[match(morn$PTID,clock$IID),9]
morn$snp2      <- clock[match(morn$PTID,clock$IID),10]

m7  <- lme(lthick ~ snp1, random = ~1|FID, data = morn, na.action = na.omit)
m8  <- lme(rthick ~ snp1, random = ~1|FID, data = morn, na.action = na.omit)
m9  <- lme(lthick ~ snp2, random = ~1|FID, data = morn, na.action = na.omit)
m10 <- lme(rthick ~ snp2, random = ~1|FID, data = morn, na.action = na.omit)

#### correlations between clock SNPs and sx
# SNPs from GWAS: snp1 = rs11121022, snp2 = rs35833281
# 2 minor alleles = 2, hets = 1, majors = 0
sx        <- read.table("swedish_sx.txt",header=T)
morn$SANS <- sx[match(morn$PTID,sx$IID),3]
morn$SAPS <- sx[match(morn$PTID,sx$IID),4]

m11 <- lme(SANS ~ snp1, random = ~1|FID, data = morn, na.action = na.omit)
m12 <- lme(SAPS ~ snp1, random = ~1|FID, data = morn, na.action = na.omit)
m13 <- lme(SANS ~ snp2, random = ~1|FID, data = morn, na.action = na.omit)
m14 <- lme(SAPS ~ snp2, random = ~1|FID, data = morn, na.action = na.omit)

# in patients
patients <- subset(morn,dx==1 | dx==3)
m15 <- lme(SANS ~ snp1, random = ~1|FID, data = patients, na.action = na.omit)
m16 <- lme(SAPS ~ snp1, random = ~1|FID, data = patients, na.action = na.omit)
m17 <- lme(SANS ~ snp2, random = ~1|FID, data = patients, na.action = na.omit)
m18 <- lme(SAPS ~ snp2, random = ~1|FID, data = patients, na.action = na.omit)

# in SZ
SZ <- subset(morn, dx==1)
m19 <- lme(SANS ~ snp1, random = ~1|FID, data = SZ, na.action = na.omit)
m20 <- lme(SAPS ~ snp1, random = ~1|FID, data = SZ, na.action = na.omit)
m21 <- lme(SANS ~ snp2, random = ~1|FID, data = SZ, na.action = na.omit)
m22 <- lme(SAPS ~ snp2, random = ~1|FID, data = SZ, na.action = na.omit)

#### correlations between clock genes and morningness
# 15 genes from GWAS = RGS16, VIP, PER2, HCRTR2 (OX2R), RASD1, 
# PER3 (VAMP3), FBXL3 (CLN5), PLCL1, APH1A (CA14), FBXL13 (FAM185A), NOL4,
# TOX4, AK5, DLX5 (SHFM1), ALG10B

source("https://bioconductor.org/biocLite.R")
biocLite("illuminaHumanv4.db")
library('illuminaHumanv4.db')

load("../RNA_GWAS/swedenclean.rdata")
clocks <- c("RGS16","VIP","PER2","HCRTR2","RASD1","PER3","FBXL3",
            "PLCL1","APH1A","FBXL13","NOL4","TOX4","AK5","DLX5","ALG10B")

xx        <- as.data.frame(illuminaHumanv4ALIAS2PROBE) # all the genes
clocks.xx <- xx[xx$alias_symbol %in% clocks,] # gwas genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% clocks.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis",probes)

df <- swedenclean[,names(swedenclean) %in% vars]
df <- merge(df,morn[,c(1,13:15,18,19)],by.x="StudyID",by.y="PTID")

RNA.morn <- matrix(nrow=length(probes),ncol=4)
for (i in 1:length(probes)){
  f <- formula(paste("morn.score", probes[i],sep="~"))
  RNA.m1 <- lme(fixed=f, random = ~1|Family, data = df, na.action = na.omit)
  tval   <- summary(RNA.m1)$tTable[2,4]
  RNA.morn[i,1] <- probes[i]
  RNA.morn[i,2] <- tval
  RNA.morn[i,3] <- sqrt(tval^2/(tval^2+summary(RNA.m1)$tTable[2,3]))
  RNA.morn[i,4] <- summary(RNA.m1)$tTable[2,5]
}

RNA.morn        <- as.data.frame(RNA.morn)
names(RNA.morn) <- c("item","t","corr","p")

#### correlations between RNA and cortical thickness
RNA.lthick <- matrix(nrow=length(probes),ncol=4)
for (i in 1:length(probes)){
  f <- formula(paste("lthick ~", probes[i],"+ Age + Sex + Diagnosis", sep=""))
  RNA.m2 <- lme(fixed=f, random = ~1|Family, data = df, na.action = na.omit)
  tval   <- summary(RNA.m2)$tTable[2,4]
  RNA.lthick[i,1] <- probes[i]
  RNA.lthick[i,2] <- tval
  RNA.lthick[i,3] <- sqrt(tval^2/(tval^2+summary(RNA.m2)$tTable[2,3]))
  RNA.lthick[i,4] <- summary(RNA.m2)$tTable[2,5]
}

RNA.lthick        <- as.data.frame(RNA.lthick)
names(RNA.lthick) <- c("item","t","corr","p")
RNA.lthick$gene   <- clocks.xx[match(RNA.lthick$item,clocks.xx$probe_id),2]

RNA.rthick <- matrix(nrow=length(probes),ncol=4)
for (i in 1:length(probes)){
  f <- formula(paste("rthick ~", probes[i],"+ Age + Sex + Diagnosis", sep=""))
  RNA.m3 <- lme(fixed=f, random = ~1|Family, data = df, na.action = na.omit)
  tval   <- summary(RNA.m3)$tTable[2,4]
  RNA.rthick[i,1] <- probes[i]
  RNA.rthick[i,2] <- tval
  RNA.rthick[i,3] <- sqrt(tval^2/(tval^2+summary(RNA.m3)$tTable[2,3]))
  RNA.rthick[i,4] <- summary(RNA.m3)$tTable[2,5]
}

RNA.rthick        <- as.data.frame(RNA.rthick)
names(RNA.rthick) <- c("item","t","corr","p")
RNA.rthick$gene   <- clocks.xx[match(RNA.rthick$item,clocks.xx$probe_id),2]

beh1 <- lme(morn.score ~ SAPS, random = ~1|FID, data = morn, na.action = na.omit)

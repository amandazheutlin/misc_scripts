# structural ROI data (freesurfer)
# swedish twins
# test against gene expression

# load data, libaries
library(nlme); library(dplyr)

workDir <- "/data/swe_gwas/ABZ/MRI_ROI_swedes/"
setwd(workDir)

# genes
load("/data/swe_gwas/ABZ/RNA_GWAS/swedenclean.rdata")
immune <- readLines("immunegenes.txt") 
immgene <- immune[immune %in% names(swedenclean)]
genes <- swedenclean[,immgene] 
genes$IID <- swedenclean[,1]

# demo
swe <- read.table("/data/ML_genetics/Schizophrenia-Zheutlin/dfs/GWASinput2-swe.txt",header=T)

# structural MRI ROI
fs_files <- list.files(pattern="FS.*\\.txt")
ROI <- lapply(fs_files[2:5], function(x) read.table(x,header=T))
for (i in 1:length(ROI)){ assign(paste0("ROI",i),as.data.frame(ROI[i])) }

######################################### overall volume against immune genes
ROI3$IID <- substr(ROI3$scanid,1,10)
sweROI <- merge(genes,ROI3,by="IID") %>% merge(swe[,1:7],by="IID")

# 25 genes against 5 phenotypes
brainvol <- c("TotalGrayVol", "lhCortexVol", "rhCortexVol", 
              "lhThickness", "rhThickness")

# total gray volume
vol_gene1 <- matrix(,nrow=length(immgene),ncol=4)
for (i in 1:length(immgene)){
  f <- formula(paste(brainvol[1],"~", immgene[i], "+ age + sex",sep=""))
  m1 <- lme(fixed=f, random = ~1|FID, data=sweROI, na.action=na.omit)
  tval <- summary(m1)$tTable[2,4]
  vol_gene1[i,1] <- immgene[i]
  vol_gene1[i,2] <- tval
  vol_gene1[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  vol_gene1[i,4] <- summary(m1)$tTable[2,5]
}

vol_gene1 <- as.data.frame(vol_gene1)
names(vol_gene1) <- c("gene","t1","corr1","p1")

# left cortex volume
vol_gene2 <- matrix(,nrow=length(immgene),ncol=4)
for (i in 1:length(immgene)){
  f <- formula(paste(brainvol[2],"~", immgene[i], "+ age + sex",sep=""))
  m1 <- lme(fixed=f, random = ~1|FID, data=sweROI, na.action=na.omit)
  tval <- summary(m1)$tTable[2,4]
  vol_gene2[i,1] <- immgene[i]
  vol_gene2[i,2] <- tval
  vol_gene2[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  vol_gene2[i,4] <- summary(m1)$tTable[2,5]
}

vol_gene2 <- as.data.frame(vol_gene2)
names(vol_gene2) <- c("gene","t2","corr2","p2")

# right cortex volume
vol_gene3 <- matrix(,nrow=length(immgene),ncol=4)
for (i in 1:length(immgene)){
  f <- formula(paste(brainvol[3],"~", immgene[i], "+ age + sex",sep=""))
  m1 <- lme(fixed=f, random = ~1|FID, data=sweROI, na.action=na.omit)
  tval <- summary(m1)$tTable[2,4]
  vol_gene3[i,1] <- immgene[i]
  vol_gene3[i,2] <- tval
  vol_gene3[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  vol_gene3[i,4] <- summary(m1)$tTable[2,5]
}

vol_gene3 <- as.data.frame(vol_gene3)
names(vol_gene3) <- c("gene","t3","corr3","p3")

# left thickness
vol_gene5 <- matrix(,nrow=length(immgene),ncol=4)
for (i in 1:length(immgene)){
  f <- formula(paste(brainvol[4],"~", immgene[i], "+ age + sex",sep=""))
  m1 <- lme(fixed=f, random = ~1|FID, data=sweROI, na.action=na.omit)
  tval <- summary(m1)$tTable[2,4]
  vol_gene4[i,1] <- immgene[i]
  vol_gene4[i,2] <- tval
  vol_gene4[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  vol_gene4[i,4] <- summary(m1)$tTable[2,5]
}

vol_gene4 <- as.data.frame(vol_gene4)
names(vol_gene4) <- c("gene","t4","corr4","p4")

# right thickness
vol_gene5 <- matrix(,nrow=length(immgene),ncol=4)
for (i in 1:length(immgene)){
  f <- formula(paste(brainvol[5],"~", immgene[i], "+ age + sex",sep=""))
  m1 <- lme(fixed=f, random = ~1|FID, data=sweROI, na.action=na.omit)
  tval <- summary(m1)$tTable[2,4]
  vol_gene5[i,1] <- immgene[i]
  vol_gene5[i,2] <- tval
  vol_gene5[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$tTable[2,3]))
  vol_gene5[i,4] <- summary(m1)$tTable[2,5]
}

vol_gene5 <- as.data.frame(vol_gene5)
names(vol_gene5) <- c("gene","t5","corr5","p5")

vol_gene <- cbind(vol_gene1,vol_gene2[2:4],vol_gene3[2:4],
                  vol_gene4[2:4],vol_gene5[2:4])

write.table(vol_gene,"swe_brainvol_immunegenes.txt",col.names=T,
            row.names=F,quote=F,sep="\t")

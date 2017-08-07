# Effect of case/control and smoking on gene expression
# Finnish twins
# 17 probes -> 7 genes

# m1: patients vs. CT/controls on gene expression
#     covar = age, sex
# m2: patients vs. CT/controls on gene expression
#     covar = age, sex, smoking status (current, former, never)
# m3: smokers vs. non-smokers on gene expression
#     covar = age, sex
# m4: smokers vs. non-smokers on gene expression
#     covar = age, sex, case/control
# m5: smoking SZ vs. non-smoking SZ vs. 
#     smoking HV vs. non-smoking HV on gene expression
#     covar = none



#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NSP_Finns/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2", "nlme", "reshape2")
invisible(lapply(libs, require, character.only = TRUE))

#### load data
load("/data/swe_gwas/ABZ/NDE1/finlandclean.rdata")
colnames(finlandclean)[2:18567] = colnames(finlandclean)[1:18566]
finlandclean     <- finlandclean[,2:18567]
finlandclean$Sex <- as.factor(finlandclean$Sex)

# lots of missing values in age + FID
# need to update using random other files
demo     <- read.table("finn-twinID-sex-age.txt",header=T,stringsAsFactors = F)
demo$age <- as.numeric(demo$age)
finlandclean$age2 <- demo[match(finlandclean$Persnb,demo$pid),6]

# USE FID NOT PAIRNB
finlandclean$FID <- demo[match(finlandclean$Persnb, demo$pid),3] %>% as.character
finlandclean[1:2,"FID"] <- finlandclean[1:2,"Pairnb"]

# add in smoking status
# (fun times with ID conversion)
# smoking file has various IID conventions and mostly missing FID
batchIDs  <- read.table("batch-conversion.txt",header=T,stringsAsFactors = F)
smoke     <- read.table("finns_smoking.txt",header=T, stringsAsFactors = F)

# split by IID convention for matching
smoke.altID <- subset(smoke, BATCH == 5 | BATCH == 6 | BATCH == 3) # altIID
smoke.corr  <- subset(smoke, BATCH == 4 | BATCH == 2) # no adjustment needed
smoke.IID   <- subset(smoke, BATCH == 1) # has IID but no FID

smoke.altID$IID <- batchIDs[match(smoke.altID$PERSON,batchIDs$IID_alt),2]
smoke.altID$FID <- batchIDs[match(smoke.altID$PERSON,batchIDs$IID_alt),1]
smoke.IID$FID   <- batchIDs[match(smoke.IID$PERSON,batchIDs$IID),1]

smoke.altID <- smoke.altID[,c("IID","FID","SMOKE_RECE","BATCH")]
smoke.IID   <- smoke.IID[,c("PERSON","FID","SMOKE_RECE","BATCH")]
smoke.corr  <- smoke.corr[,c("PERSON","FAMILY","SMOKE_RECE","BATCH")]

names.cols <- c("IID", "FID", "smoke", "batch")
names(smoke.altID) <- names.cols
names(smoke.IID)   <- names.cols
names(smoke.corr)  <- names.cols

smoke.use <- rbind(smoke.altID, smoke.corr,smoke.IID)

finlandclean$smoke <- smoke.use[match(finlandclean$Persnb, smoke.use$IID),3] 
finlandclean$smoke <- ifelse(finlandclean$smoke==9,NA,finlandclean$smoke) # 9 = missing
           # missing 19 of 73

# NSP genes
genes        <- c("ILMN_1658472","ILMN_2398388","ILMN_1862217","ILMN_1767816",
                  "ILMN_2320349","ILMN_1797804","ILMN_1653728","ILMN_1735180",
                  "ILMN_1731788","ILMN_1737252","ILMN_1657607","ILMN_1811954",
                  "ILMN_2298888","ILMN_1774945","ILMN_1809193","ILMN_2404512",
                  "ILMN_1714417")
vars         <- c("Persnb","FID","Status","age2", "Sex", "smoke", genes)
df           <- finlandclean[,colnames(finlandclean) %in% vars]
probes       <- names(df)[4:15]

############ stats - mixed anovas
# m1: DV = gene expression
# IV = case vs. control+CT
# covar = age, sex

m1 <- lapply(probes, function(x) {
  lme(eval(substitute(i ~ Status + age2 + Sex, list(i = as.name(x)))), 
      random = ~1|FID, data = df, na.action = na.omit)
})

# stats
m1_stats <- lapply(m1, function(x) summary(x)$tTable)

# naming genes
getResponse <- function(formula) {
  tt <- summary(formula)$terms
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response] 
}

# list of overall effects
NSPm1 = NULL
for (i in 1:12) {
  temp <- as.data.frame(m1_stats[[i]])
  temp$marker <- getResponse(m1[[i]])
  NSPm1 <- rbind(NSPm1,temp[2,])
}

colnames(NSPm1)[4:5] <- c("tvalue","pvalue")
NSPm1                <- NSPm1[order(NSPm1$pvalue),]

write.table(NSPm1,"cc_age-sex.NSPm1.txt",sep="\t",
            col.names=T,row.names=F,quote=F)

# m2: DV = gene expression
# IV = case vs. control+CT
# covar = age, sex, smoked recently (yes/no)

m2 <- lapply(probes, function(x) {
  lme(eval(substitute(i ~ Status + age2 + Sex + smoke, list(i = as.name(x)))), 
      random = ~1|FID, data = df, na.action = na.omit)
})

# stats
m2_stats <- lapply(m2, function(x) summary(x)$tTable)

# list of overall effects
NSPm2 = NULL
for (i in 1:12) {
  temp <- as.data.frame(m2_stats[[i]])
  temp$marker <- getResponse(m2[[i]])
  NSPm2 <- rbind(NSPm2,temp[2,])
}

colnames(NSPm2)[4:5] <- c("tvalue","pvalue")
NSPm2                <- NSPm2[order(NSPm2$pvalue),]

write.table(NSPm2,"cc_age-sex-smoke.NSPm2.txt",sep="\t",
            col.names=T,row.names=F,quote=F)


# m3: DV = gene expression
# IV = smokers vs. non-smokers
# covar = age, sex

m3 <- lapply(probes, function(x) {
  lme(eval(substitute(i ~ smoke + age2 + Sex, list(i = as.name(x)))), 
      random = ~1|FID, data = df, na.action = na.omit)
})

# stats
m3_stats <- lapply(m3, function(x) summary(x)$tTable)

# list of overall effects
NSPm3 = NULL
for (i in 1:12) {
  temp <- as.data.frame(m3_stats[[i]])
  temp$marker <- getResponse(m3[[i]])
  NSPm3 <- rbind(NSPm3,temp[2,])
}

colnames(NSPm3)[4:5] <- c("tvalue","pvalue")
NSPm3                <- NSPm3[order(NSPm3$pvalue),]

write.table(NSPm3,"smoke_age-sex.NSPm3.txt",sep="\t",
            col.names=T,row.names=F,quote=F)


# m4: DV = gene expression
# IV = smokers vs. non-smokers
# covar = age, sex, case/control

m4 <- lapply(probes, function(x) {
  lme(eval(substitute(i ~ smoke + age2 + Sex + Status, list(i = as.name(x)))), 
      random = ~1|FID, data = df, na.action = na.omit)
})

# stats
m4_stats <- lapply(m4, function(x) summary(x)$tTable)

# list of overall effects
NSPm4 = NULL
for (i in 1:12) {
  temp <- as.data.frame(m4_stats[[i]])
  temp$marker <- getResponse(m4[[i]])
  NSPm4 <- rbind(NSPm4,temp[2,])
}

colnames(NSPm4)[4:5] <- c("tvalue","pvalue")
NSPm4                <- NSPm4[order(NSPm4$pvalue),]

write.table(NSPm4,"smoke_age-sex-cc.NSPm4.txt",sep="\t",
            col.names=T,row.names=F,quote=F)


# m5: DV = gene expression
# IV = SZ smokers vs. SZ non-smokers vs. 
#      HV smokers vs. HV non-smokers
# covar = none

df2 <- df %>% 
          filter(complete.cases(smoke)) %>%
              group_by(Status,smoke) %>%
              mutate(.,smoke_cc = paste0(Status,smoke))

df2$smoke_cc         <- as.factor(df2$smoke_cc)
levels(df2$smoke_cc) <- c('HVns', 'HVs', 'SZns', 'SZs')

m5 <- lapply(probes, function(x) {
  lme(eval(substitute(i ~ smoke_cc, list(i = as.name(x)))), 
      random = ~1|FID, data = df2, na.action = na.omit)
})

# stats
m5_stats <- lapply(m5, function(x) summary(x)$tTable)

# list of overall effects
NSPm5 = NULL
for (i in 1:12) {
  temp <- as.data.frame(m5_stats[[i]])
  temp$marker <- getResponse(m5[[i]])
  NSPm5 <- rbind(NSPm5,temp[2,])
}

colnames(NSPm5)[4:5] <- c("tvalue","pvalue")
NSPm5                <- NSPm5[order(NSPm5$pvalue),]

write.table(NSPm5,"smokecc.NSPm5.txt",sep="\t",
            col.names=T,row.names=F,quote=F)





################################################
# checking N for twin batch
df$twin.no <- demo[match(df$Persnb,demo$pid),1]
df$batch   <- substring(df$twin.no,1,3)
table(df$batch) # 12, 26, 8, 8, 17

df$batch2 <- batchIDs[match(df$Persnb,batchIDs$IID),5]
table(df$batch2) # 12, 25, 10, 8, 17




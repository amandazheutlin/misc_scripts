# check 2 SNPs in NDE1 against genome-wide gene expression
# Finnish twins 

# load data, libraries
libs <- c("nlme","dplyr")
invisible(lapply(libs, require, character.only = TRUE))

workDir <- "/data/swe_gwas/ABZ/NDE1/"
setwd(workDir)

load("finlandclean.rdata")
colnames(finlandclean)[2:18567] <- colnames(finlandclean)[1:18566]
finlandclean <- finlandclean[,2:18567]

demo     <- read.table("demo_finns.txt",header=T,stringsAsFactors = F)
demo$age <- as.numeric(demo$age)
finlandclean$age2 <- demo[match(finlandclean$Persnb,demo$pid),6]

geno        <- read.table("NDE1.ped",header=F,sep="\t")
geno        <- geno[,c(1:2,7:8)]
names(geno) <- c("FID","IID","rs2242549","rs1050162")

# update finlandclean FIDs that are missing!
# USE FID NOT PAIRNB
finlandclean$FID <- geno[match(finlandclean$Persnb, geno$IID),1] %>% as.character

# snps = rs2242549 (snp1), rs1050162 (snp2)
# snps were coded 1 1, 1 2, 2 2; 1 = minor; 2 = major
# minor allele = G (snp1), C (snp2); major allele = T (both)
# additive = 1, 2, 3
# recessive = 1, 2, 2 (name = rs rec)
# dominant = 1, 1, 2 (name = rs dom)
NDE1 <- finlandclean
NDE1$rs2242549 <- geno[match(NDE1$Persnb,geno$IID),3] 
NDE1$rs2242549rec <- NDE1$rs2242549
NDE1$rs2242549dom <- NDE1$rs2242549
NDE1$rs1050162 <- geno[match(NDE1$Persnb,geno$IID),4] 
NDE1$rs1050162rec <- NDE1$rs1050162
NDE1$rs1050162dom <- NDE1$rs1050162

NDE1$rs2242549 <- ifelse(NDE1$rs2242549=="1 1", 1, ifelse(NDE1$rs2242549=="1 2",2,3)) # N(1,2,3) = 6,38,26
NDE1$rs2242549rec <- ifelse(NDE1$rs2242549rec=="1 1", 1, 2) # N(1,2) = 6,64
NDE1$rs2242549dom <- ifelse(NDE1$rs2242549dom=="2 2", 2, 1) # N(1,2) = 44,26
NDE1$rs1050162 <- ifelse(NDE1$rs1050162=="1 1", 1, ifelse(NDE1$rs1050162=="1 2",2,3)) # N(1,2,3) = 10,37,23
NDE1$rs1050162rec <- ifelse(NDE1$rs1050162rec=="1 1", 1, 2) # N(1,2) = 10,60
NDE1$rs1050162dom <- ifelse(NDE1$rs1050162dom=="2 2", 2, 1) # N(1,2) = 47,23

genes <- names(NDE1)[8:18566] # 18,559 markers

####################################### residual values for all genes
# effects of sex, diagnosis ('Analysis'), and age (use 'age2')
# Analysis = case/co-twin/control

df        <- NDE1[complete.cases(NDE1$age2),]
df        <- df[complete.cases(df$rs2242549),]
residuals <- lapply(genes, function(x) {
  lm(eval(substitute(i ~ Sex + Analysis + age2, list(i = as.name(x)))), 
      data = df, na.action = na.omit)
})

res_vals     <- lapply(residuals, function(x) residuals(x)) %>% as.data.frame()
resdf        <- cbind(df[,c(1,18568,18569,18571,18572,18574)],res_vals)
names(resdf) <- c("IID","FID","rs2242549","rs2242549dom",
                  "rs1050162","rs1050162dom",genes)

####################################### snps x expression (residuals)
# effect of snp on all genes (N = 18,559 genes)
# residual values excluding effects of sex, diagnosis ('Analysis')
# random effect = family 
# using lmeControl to deal with convergence issues

#########################################################
# snp1 - additive model
models_snp1 <- lapply(genes, function(x) {
  lme(eval(substitute(i ~ rs2242549, list(i = as.name(x)))), 
      random = ~1|FID, data = resdf, na.action = na.omit,
      control=lmeControl(returnObject=T, opt="optim"))
})

# stats
model_stats_snp1 <- lapply(models_snp1, function(x) summary(x)$tTable)

# naming genes
getResponse <- function(formula) {
  tt <- summary(formula)$terms
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response] 
}

# list of overall effects
snp1 = NULL
for (i in 1:18559) {
  temp <- as.data.frame(model_stats_snp1[[i]])
  temp$marker <- getResponse(models_snp1[[i]])
  snp1 <- rbind(snp1,temp[2,])
}

colnames(snp1)[4:5] <- c("tvalue","pvalue")
snp1 <- snp1[order(snp1$pvalue),]

write.table(snp1,"RNA-GWAS-snp1-add-res-agecov-june16.txt",
            sep="\t",col.names=T,row.names=F,quote=F)

rm(models_snp1)
rm(model_stats_snp1)
rm(snp1)
#########################################################
# snp1 - dominant model
models_snp1dom <- lapply(genes, function(x) {
  lme(eval(substitute(i ~ rs2242549dom, list(i = as.name(x)))), 
      random= ~1|FID, data = resdf, na.action = na.omit,
      control=lmeControl(returnObject=T, opt="optim"))
})

# stats
model_stats_snp1dom <- lapply(models_snp1dom, function(x) summary(x)$tTable)

# list of overall effects
snp1_dom = NULL
for (i in 1:18559) {
  temp <- as.data.frame(model_stats_snp1dom[[i]])
  temp$marker <- getResponse(models_snp1dom[[i]])
  snp1_dom <- rbind(snp1_dom,temp[2,])
}

colnames(snp1_dom)[4:5] <- c("tvalue","pvalue")
snp1_dom <- snp1_dom[order(snp1_dom$pvalue),]

write.table(snp1_dom,"RNA-GWAS-snp1-dom-res-agecov-june16.txt",
            sep="\t",col.names=T,row.names=F,quote=F)

rm(models_snp1dom)
rm(model_stats_snp1dom)
rm(snp1_dom)
#########################################################
# snp2 - additive model
models_snp2 <- lapply(genes, function(x) {
  lme(eval(substitute(i ~ rs1050162, list(i = as.name(x)))), 
      random= ~1|FID, data = resdf, na.action = na.omit,
      control=lmeControl(returnObject=T, opt="optim"))
})

# stats
model_stats_snp2 <- lapply(models_snp2, function(x) summary(x)$tTable)

# list of overall effects
snp2 = NULL
for (i in 1:18559) {
  temp <- as.data.frame(model_stats_snp2[[i]])
  temp$marker <- getResponse(models_snp2[[i]])
  snp2 <- rbind(snp2,temp[2,])
}

names(snp2)[4:5] <- c("tvalue","pvalue")
snp2 <- snp2[order(snp2$pvalue),]

write.table(snp2,"RNA-GWAS-snp2-add-res-agecov-june16.txt",
            sep="\t",col.names=T,row.names=F,quote=F)

rm(models_snp2)
rm(model_stats_snp2)
rm(snp2)
#########################################################
# snp2 - dominant model
models_snp2dom <- lapply(genes, function(x) {
  lme(eval(substitute(i ~ rs1050162dom, list(i = as.name(x)))), 
      random= ~1|FID, data = resdf, na.action = na.omit,
      control=lmeControl(returnObject=T, opt="optim"))
})

# stats
model_stats_snp2dom <- lapply(models_snp2dom, function(x) summary(x)$tTable)

# list of overall effects
snp2_dom = NULL
for (i in 1:18559) {
  temp <- as.data.frame(model_stats_snp2dom[[i]])
  temp$marker <- getResponse(models_snp2dom[[i]])
  snp2_dom <- rbind(snp2_dom,temp[2,])
}

names(snp2_dom)[4:5] <- c("tvalue","pvalue")
snp2_dom <- snp2_dom[order(snp2_dom$pvalue),]

write.table(snp2_dom,"RNA-GWAS-snp2-dom-res-agecov-june16.txt",
            sep="\t",col.names=T,row.names=F,quote=F)

rm(models_snp2dom)
rm(model_stats_snp2dom)
rm(snp2_dom)

###########################################################
# tvab
famID <- unique(df$FID)
twin1 <- df[!duplicated(df$FID),]
twin2 <- df[!(df$Persnb %in% twin1$Persnb),]

###########################################################
###########################################################
####################################### snps x expression
# effect of snp on all genes (N = 18,559 genes)
# covar = sex, diagnosis ('Analysis')
# random effect = family 
# Analysis = case/co-twin/control

# #########################################################
# # snp1 - additive model
# models_snp1 <- lapply(genes, function(x) {
#   lme(eval(substitute(i ~ rs2242549 + Sex + Analysis, list(i = as.name(x)))), 
#       random = ~1|Pairnb, data = NDE1, na.action = na.omit)
# })
# 
# # stats
# model_stats_snp1 <- lapply(models_snp1, function(x) summary(x)$tTable)
# 
# # naming genes
# getResponse <- function(formula) {
#   tt <- summary(formula)$terms
#   vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
#   response <- attr(tt, "response") # index of response var
#   vars[response] 
# }
# 
# # list of overall effects
# snp1 = NULL
# for (i in 1:18559) {
#   temp <- as.data.frame(model_stats_snp1[[i]])
#   temp$marker <- getResponse(models_snp1[[i]])
#   snp1 <- rbind(snp1,temp[2,])
# }
# 
# colnames(snp1)[4:5] <- c("tvalue","pvalue")
# snp1 <- snp1[order(snp1$pvalue),]
# 
# write.table(snp1,"RNA-GWAS-finns-snp1-add.txt",
#             sep="\t",col.names=T,row.names=F,quote=F)
# 
# #########################################################
# # snp1 - recessive model
# models_snp1rec <- lapply(genes, function(x) {
#   lme(eval(substitute(i ~ rs2242549rec + Sex + Analysis, list(i = as.name(x)))), 
#       random= ~1|Pairnb, data = NDE1, na.action = na.omit)
# })
# 
# # stats
# model_stats_snp1rec <- lapply(models_snp1rec, function(x) summary(x)$tTable)
# 
# # list of overall effects
# snp1_rec = NULL
# for (i in 1:18559) {
#   temp <- as.data.frame(model_stats_snp1rec[[i]])
#   temp$marker <- getResponse(models_snp1rec[[i]])
#   snp1_rec <- rbind(snp1_rec,temp[2,])
# }
# 
# colnames(snp1_rec)[4:5] <- c("tvalue","pvalue")
# snp1_rec <- snp1_rec[order(snp1_rec$pvalue),]
# 
# write.table(snp1_rec,"RNA-GWAS-finns-snp1-rec.txt",
#             sep="\t",col.names=T,row.names=F,quote=F)
# 
# #########################################################
# # snp2 - additive model
# models_snp2 <- lapply(genes, function(x) {
#   lme(eval(substitute(i ~ rs1050162 + Sex + Analysis, list(i = as.name(x)))), 
#       random= ~1|Pairnb, data = NDE1, na.action = na.omit)
# })
# 
# # stats
# model_stats_snp2 <- lapply(models_snp2, function(x) summary(x)$tTable)
# 
# # list of overall effects
# snp2 = NULL
# for (i in 1:18559) {
#   temp <- as.data.frame(model_stats_snp2[[i]])
#   temp$marker <- getResponse(models_snp2[[i]])
#   snp2 <- rbind(snp2,temp[2,])
# }
# 
# names(snp2)[4:5] <- c("tvalue","pvalue")
# snp2 <- snp2[order(snp2$pvalue),]
# 
# write.table(snp2,"RNA-GWAS-finns-snp2-add.txt",
#             sep="\t",col.names=T,row.names=F,quote=F)
# 
# #########################################################
# # snp2 - recessive model
# models_snp2rec <- lapply(genes, function(x) {
#   lme(eval(substitute(i ~ rs1050162rec + Sex + Analysis, list(i = as.name(x)))), 
#       random= ~1|Pairnb, data = NDE1, na.action = na.omit)
# })
# 
# # stats
# model_stats_snp2rec <- lapply(models_snp2rec, function(x) summary(x)$tTable)
# 
# # list of overall effects
# snp2_rec = NULL
# for (i in 1:18559) {
#   temp <- as.data.frame(model_stats_snp2rec[[i]])
#   temp$marker <- getResponse(models_snp2rec[[i]])
#   snp2_rec <- rbind(snp2_rec,temp[2,])
# }
# 
# names(snp2_rec)[4:5] <- c("tvalue","pvalue")
# snp2_rec <- snp2_rec[order(snp2_rec$pvalue),]
# 
# write.table(snp2_rec,"RNA-GWAS-finns-snp2-rec.txt",
#             sep="\t",col.names=T,row.names=F,quote=F)
# 

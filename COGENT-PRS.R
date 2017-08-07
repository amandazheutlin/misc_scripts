# PsychChip data from NAPLS - SNPs
# COGENT PRS scores reflecting genetic variants associated with cognition
# correlations with risk, conversion, SZ PRS, neuropsych

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/COGENT/"
setwd(workdir)

libs <- c("dplyr", "ggplot2", "data.table", "gdata")
invisible(lapply(libs, require, character.only = TRUE))

######## load in data
# file with demo, PRS, ancestry PCs
# generated with genos_conversion.R (in NAPLS/PSYCH-CHIP/)
demo.prs <- read.table("../PSYCH-CHIP/napls-SZrisk-PCs-conv.txt", 
                       header=T, stringsAsFactors = F, sep="\t")

# file with COGENT summary data
# it is a bit cumbersome (.66 GB)
# cog.sum     <- fread("cogent.hrc.meta.chr.bp.rsid.assoc.full", header=T, stringsAsFactors = F)
# cog.input   <- cog.sum[!(cog.sum$RSID == "-"),] %>%       # get rid of nonsense rows
#                dplyr::select(RSID, EFFECT_ALLELE, BETA) 
# cog.input.5 <- cog.sum[!(cog.sum$RSID == "-"),] %>%       # get rid of nonsense rows
#                dplyr::select(RSID, EFFECT_ALLELE, BETA, PVAL) 
# cog.input.5 <- cog.input.5[cog.input.5$PVAL<.5,] %>%       
#                dplyr::select(RSID, EFFECT_ALLELE, BETA)
# 
# write.table(cog.input, "allsnps-COGprs.txt", sep="\t", col.names=F, row.names=F, quote=F)
# write.table(cog.input.5, "pT5-COGprs.txt", sep="\t", col.names=F, row.names=F, quote=F)
# 
# rm(cog.sum)
# rm(cog.input)
# rm(cog.input.5)
# 
# generate scores using input file in plink
# run in PSYCH-CHIP/dosage/
# plink --dosage /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/dosagelist.txt list format=2 --fam dos_scz_scz1_mix_la-qc.hg19.ch.fl.chr9_141_144.out.dosage.fam --score /data/swe_gwas/ABZ/NAPLS/COGENT/allsnps-COGprs.txt header sum --out /data/swe_gwas/ABZ/NAPLS/COGENT/napls_COGprs_allsnps
# plink --dosage /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/dosagelist.txt list format=2 --fam dos_scz_scz1_mix_la-qc.hg19.ch.fl.chr9_141_144.out.dosage.fam --score /data/swe_gwas/ABZ/NAPLS/COGENT/pT5-COGprs.txt header sum --out /data/swe_gwas/ABZ/NAPLS/COGENT/napls_COGprs_pT5

# file with COGENT scores
cog.prs <- read.table("napls_COGprs_allsnps.profile", header=T, stringsAsFactors = F)
cog.pT5 <- read.table("napls_COGprs_pT5.profile", header=T, stringsAsFactors = F)
demo.prs$cogprs <- cog.prs[match(demo.prs$IID, cog.prs$IID), "SCORESUM"]
demo.prs$cogpT5 <- cog.pT5[match(demo.prs$IID, cog.pT5$IID), "SCORESUM"]

# file with neuropsych data
# needs ID conversion
np <- read.csv("NAPLSII_imputation_and_factors_imputation5-YC.csv", header=T, stringsAsFactors = F)

id.conv             <- read.xls("../PSYCH-CHIP/NAPLS2_IDconversion_edited.xls", stringsAsFactors = F)
np$SiteSubjID       <- paste0(np$Site, np$Subject)
demo.prs$site       <- id.conv[match(demo.prs$IID, id.conv$Collaborator.Sample.ID), "NAPLS2.Site.Number"]
demo.prs$subj       <- id.conv[match(demo.prs$IID, id.conv$Collaborator.Sample.ID), "NAPLS2.Site.ID"]
demo.prs$SiteSubjID <- paste0(demo.prs$site,demo.prs$subj)

allvars <- merge(demo.prs, np, by="SiteSubjID")

# graph PCs
# ggplot(demo.prs, aes(x=PC3,y=PC2)) + 
#   geom_point(size=3) +
#   theme(axis.title = element_text(size=15, face="bold"),
#         axis.text = element_text(size=13, color="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(color = "black"), 
#         panel.grid.minor = element_line(color = "black"),
#         legend.position = "bottom",
#         legend.text = element_text(size=13, color="black"),
#         legend.title = element_blank(),
#         legend.key = element_blank(),
#         axis.line.x = element_line(color = "black"),
#         axis.line.y = element_line(color = "black")) +
#   ylab("PC2") +
#   xlab("PC3")

######## stats

# regress out ancestry from cognitive PRS
# use residuals in all analyses
cog.anc            <- lm(cogprs ~ PC1 + PC2 + PC3 + PC4 + PC5, data = demo.prs, na.action = na.omit)
cog.anc.np         <- lm(cogprs ~ PC1 + PC2 + PC3 + PC4 + PC5, data = allvars, na.action = na.omit)
pT5.anc            <- lm(cogpT5 ~ PC1 + PC2 + PC3 + PC4 + PC5, data = demo.prs, na.action = na.omit)
demo.prs$cog.resid <- residuals(cog.anc)
demo.prs$pT5.resid <- residuals(pT5.anc)
allvars$cog.resid  <- residuals(cog.anc.np)

# np factor 1 x cognition PRS
# factor 1 is exec fx / visual abilities
cog.F1 <- lm(Factor1 ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.F1)

# np correlation matric
npcor <- allvars[,c(60:77,79:174)]
cormat <- cor(npcor$cog.resid,npcor, use="pairwise.complete.obs") %>% as.data.frame %>% t
cor.df <- cbind(rownames(cormat),cormat[,1]) %>% as.data.frame
cor.df <- cor.df[order(cor.df$V2),]
cor.df[104:114,]

# cog.trail1 <- lm(hvlttrail1 ~ cog.resid + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
# summary(cog.trail1)

# WASI Vocabulary
cog.wasi.v <- lm(wasiraw ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.wasi.v)

# WASI Block Design
cog.wasi.bd <- lm(wasiblraw ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.wasi.bd)

# HVLT-R
cog.hvlt <- lm(hvlttotal ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.hvlt)

# LNS
cog.lns <- lm(lnstotal ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.lns)

# Factor 4 (declarative memory)
cog.F4 <- lm(Factor4 ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.F4)

# Factor 2 (verbal abilities)
cog.F2 <- lm(Factor2 ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.F2)

# fluency
cog.flu <- lm(fluencycor ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.flu)

# WASI IQ
cog.iq <- lm(wasiiq ~ cog.resid + dx + sex.x + age + age*sex.x + site, data = allvars, na.action=na.omit)
summary(cog.iq)

# (resid) cognition PRS x SZ PRS
# lots of covars
# turns out this doesn't really make sense and 
# no one ever correlates two PRS within individuals
# cogSZ.pTall <- lm(pTall ~ cog.resid + sex + age + age*sex + site, data = demo.prs, na.action=na.omit)
# summary(cogSZ.pTall)
# cogSZ.pT5   <- lm(pTall ~ pT5.resid + sex + age + age*sex + site, data = demo.prs, na.action=na.omit)
# summary(cogSZ.pT5)

# risk x (resid) cognition PRS
demo.prs$risk.coded <- ifelse(demo.prs$risk=="Control", 0, 1)
demo.prs$risk.coded <- factor(demo.prs$risk.coded)
cog.risk  <- glm(risk.coded ~ cog.resid + sex + age + age*sex + site, family = binomial(link="logit"), data = demo.prs, na.action=na.omit)
summary(cog.risk)

# conversion x (resid) cognition PRS in UHR
demo.prs.uhr          <- demo.prs[demo.prs$risk == "Prodromal",]
demo.prs.uhr$dx.coded <- ifelse(demo.prs.uhr$dx=="Clinical High Risk-Not Converted", 0, 1)
demo.prs.uhr$dx.coded <- factor(demo.prs.uhr$dx.coded)
cog.dx <- glm(dx.coded ~ cog.resid + sex + age + age*sex + site, family = binomial(link="logit"), data = demo.prs.uhr, na.action=na.omit)
summary(cog.dx)


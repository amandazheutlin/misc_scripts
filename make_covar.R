## make covariates csv file for ENIGMA
## swedish twins -- subject file from Yoon

######## housekeeping
workdir <- "/data/swe_gwas/ABZ/ENIGMA/"
setwd(workdir)

library(dplyr)

######## load in data
# subject list
subs <- read.table("List_subjects.txt", header=F, stringsAsFactors = F)

# correct IIDs to not include date
subs$IID <- substr(subs$V1,1,10)

# load in previously generated swe twins file
demo <- read.table("swedishNP_PCA.txt", header=T, stringsAsFactors = F)

# are subs in demo file? 
table(subs$IID %in% demo$IID) # 300/300 IDs present!

######## add columns to subject file as specified for ENIGMA
######## SCHIZOPHRENIA
# SubjID (IID), FID, age, sex (1 = male), affected (1 = SZ, 0 = all other),
# cotwin (1 = cotwin, 0 = all other), zygosity (MZ = 1, DZ = 0),
# affected_otherwise (1 = not healthy or SZ/BP, 0 = all other)
# diagnosis (standard), NO SITE INFO

enigma <- demo %>% filter(demo$IID %in% subs$IID) %>%
          dplyr::select(IID, FID, age, sex, zyg, subjtype)

names(enigma)             <- c("SubjID", "FID", "Age", "Sex", "Zygosity", "Diagnosis")
enigma$Affected           <- ifelse(enigma$Diagnosis == "sch-proband", 1, 0)
enigma$Cotwin             <- ifelse(enigma$Diagnosis == "sch-cotwin", 1, 0)
enigma$Affected_otherwise <- ifelse(enigma$Diagnosis == "dep-proband", 1, 
                                    ifelse(enigma$Diagnosis == "opsy-proband", 1, 0))
enigma$Zygosity           <- ifelse(enigma$Zygosity == "MZ", 1, 0)

# check for missing data and correct calculations
table(is.na(enigma))     # 2 values missing in df
table(is.na(enigma$Age)) # none missing
table(is.na(enigma$Sex)) # 2 subs missing (found the missing values!)
table(enigma$Zygosity)   # correct transformation
table(enigma$Affected)   # matches dx var
table(enigma$Cotwin)     # matches dx var
table(enigma$Affected_otherwise) # matches dx var

# manually update missing sex info from any twin info
enigma[is.na(enigma$Sex),] 
enigma[is.na(enigma$Sex), "Sex"] <- 9 # R doesn't like NAs

# i have a master sheet with all this info
# swedish_id_MASTER.xlsx
enigma$Sex[enigma$SubjID == 2006600857] <- 2
enigma$Sex[enigma$SubjID == 2006601157] <- 1

# write schizophrenia table
write.csv(enigma, "Covariates_SZ.csv", row.names=F, quote=F)


######## add columns to subject file as specified for ENIGMA
######## BIPOLAR
# SubjID (IID), FID, age, sex (1 = male), affected (1 = BP, 0 = all other),
# cotwin (1 = cotwin, 0 = all other), zygosity (MZ = 1, DZ = 0),
# affected_otherwise (1 = not healthy or SZ/BP, 0 = all other)
# diagnosis (standard), NO SITE INFO

enigma <- demo %>% filter(demo$IID %in% subs$IID) %>%
  dplyr::select(IID, FID, age, sex, zyg, subjtype)

names(enigma)             <- c("SubjID", "FID", "Age", "Sex", "Zygosity", "Diagnosis")
enigma$Affected           <- ifelse(enigma$Diagnosis == "bip-proband", 1, 0)
enigma$Cotwin             <- ifelse(enigma$Diagnosis == "bip-cotwin", 1, 0)
enigma$Affected_otherwise <- ifelse(enigma$Diagnosis == "dep-proband", 1, 
                                    ifelse(enigma$Diagnosis == "opsy-proband", 1, 0))
enigma$Zygosity           <- ifelse(enigma$Zygosity == "MZ", 1, 0)

# check for missing data and correct calculations
table(is.na(enigma))     # 2 values missing in df --> sex
enigma$Sex[enigma$SubjID == 2006600857] <- 2
enigma$Sex[enigma$SubjID == 2006601157] <- 1

table(enigma$Affected)   # matches dx var
table(enigma$Cotwin)     # matches dx var
table(enigma$Affected_otherwise) # matches dx var

# write bipolar table
write.csv(enigma, "Covariates_BP.csv", row.names=F, quote=F)



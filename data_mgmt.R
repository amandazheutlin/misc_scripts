# data for Cathy Zhu - Fall 2016
# neuropsych in Swedish twins

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/CZ-thesis/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2", "nlme")
invisible(lapply(libs, require, character.only = TRUE))

#### vars
# IID, FID, subjtype, pairtype, age, education
# Kessler social support, SANS, SAPS
# WAIS performance, WAIS verbal, WAIS total
# WAIS perf = block design; verbal = vocab; no total
# VR I, VR II, digit span, visual span, verbal WM, TMT
demo <- read.table("../swedish_RNAseq/swedish_demo.txt", header=T, stringsAsFactors = F)
kess <- read.csv("kessler_social_support.txt", header=T, sep="\t", stringsAsFactors = F)
SANS <- read.csv("SANS.txt", header=T, sep="\t", stringsAsFactors = F)
SAPS <- read.csv("SAPS.txt", header=T, sep="\t", stringsAsFactors = F)
NP   <- read.csv("np_all.txt", header=T, sep="\t", stringsAsFactors = F)

# merge vars
names(demo)[2] <- "PTID"
data <- demo %>% full_join(kess,by="PTID") %>% 
  full_join(SANS,by="PTID") %>% 
  full_join(SAPS, by="PTID") %>%
  full_join(NP, by="PTID") 

data.clean <- data[,c(1,2,4,6,3,7,5,94:106,15:26,31:55,60:93)]
names(data.clean)[1:20] <- c("FID","IID","pairtype","subjtype",
                            "sex","age","zygosity","yrs_edu",
                            "vocab","blck","VR-I","VR-II",
                            "digitspan","TMT-A","TMT-B","spatialspan",
                            "verbalfluency","animals","fruits","vegetables")

write.table(data.clean,"CZ_thesis_data.txt",col.names=T,
            row.names=F,quote=F,sep="\t")


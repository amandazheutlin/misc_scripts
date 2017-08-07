# demo info for RNAseq subs

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/swedish_RNAseq/"
setwd(workdir)

libs <- c("dplyr", "ggplot2", "lubridate", "xlsx")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
# demo is from swedish_id_MASTER.xlsx in genetics imputation folder (local)
# antibody-results is from Bob
demo <- read.table("swedish_demo.txt",header=T,stringsAsFactors = F)
subs <- read.table("antibody-results.txt",header=T,stringsAsFactors = F)

# calculate age (it's wrong??)
# birth month/year - testdate
demo$test.date <- format(as.Date(demo$testdate, format="%m/%d/%y"), format="%m/%d/%Y") %>% mdy()
demo$bday.date <- paste(demo$mobirth,"/01/",demo$yrbirth,sep="") %>% mdy()

demo$age.calc  <- demo$test.date - demo$bday.date 
demo$age.calc  <- demo$age.calc / 365.2425 # leapdays
demo$age.calc  <- as.numeric(demo$age.calc)

# duplicate entries
dup.ID     <- demo[duplicated(demo$pid),2]
dups       <- demo[demo$pid %in% dup.ID,] # no obvious inconsistencies
demo.clean <- demo[!duplicated(demo$pid),]

# add dx (coded A, B, C), zygosity, age, sex
# A-I: BP-CT, BP, HV, MDD-CT, MDD, opsy-CT, opsy, SZ-CT, SZ
demo.clean$dx.coded         <- as.factor(demo.clean$subjtype)
levels(demo.clean$dx.coded) <- LETTERS[1:9]

subs.demo <- merge(subs, demo.clean[,c("pairnum","pid","subjtype","dx.coded","zyg","age.calc","Sex")],
                   by.x="iid",by.y="pid")

# add birth location
birth                <- read.csv("birthplace.txt", header=T, stringsAsFactors = F, sep="\t")
subs.demo$birthplace <- birth[match(subs.demo$iid,birth$PTID),2]

# demo + results
names(subs.demo)[20:25] <- c("FID","subjtype","dx.coded","zygosity","age","sex")
write.table(subs.demo,"RNAseq_swe_demo2.txt", col.names=T,
            row.names=F, quote=F, sep="\t")
write.xlsx(subs.demo,"RNAseq_swe_demo.xlsx", col.names=T, row.names=F)

# file for bob
# jhunumber, iid, fid, age, sex, birthplace 
# zygosity, pairtype, subjtype
names(sub.demo)
subs.demo$pairtype <- demo.clean[match(subs.demo$iid,demo.clean$pid),4]

demo.wr <- subs.demo %>% dplyr::select(jhuno,iid,FID,age,sex,birthplace,zygosity,pairtype,subjtype)
write.csv(demo.wr, file = "RNAseq_swe_demo3.csv", row.names=F, quote=F)

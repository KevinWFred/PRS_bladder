#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))

phenotype=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
dups=read.csv("../data/DUPS_PLCOdata.csv")
table(dups$plco_id %in% phenotype$study_pid)
for (i in 1:ncol(phenotype))
{
  idx=which(phenotype[,i]=="")
  if (length(idx)>0) phenotype[idx,i]=NA
  idx=which(phenotype[,i]=="MISSING")
  if (length(idx)>0) phenotype[idx,i]=NA
  if (any(dups$plco_id %in% phenotype[,i]))
  {
    print(table(dups$plco_id %in% phenotype[,i]))
  }
}
bladder24pc=read.csv("../data/bladder_individual_24PCs.csv")
table(bladder24pc$study_pid %in% phenotype$study_pid,useNA = "ifany")
# FALSE  TRUE 
# 214 20275
table(bladder24pc$pid %in% phenotype$pid,useNA = "ifany")
# FALSE  TRUE 
# 41 20448
table(bladder24pc$gwas_id %in% phenotype$gwas_id,useNA = "ifany")
# FALSE  TRUE 
# 41 20448

#library("sas7bdat")
#dat=read.sas7bdat("../data/package-plco-957/Bladder File/plco_957_033022.sas7bdat")
#write.csv(dat,file="../data/package-plco-957/Bladder File/plco_957_033022.csv",row.names = F)
dat=as.data.frame(fread("../data/package-plco-957/Bladder File/plco_957_033022.csv"))
table(dups$plco_id %in% dat$plco_id)
idx=match(dups$plco_id,dat$plco_id)
table(dat$exc_in_tgwas_population[idx],useNA = "ifany")
# 1 
# 68228
table(dups$in_prev_analysis,useNA = "ifany")
# 0     1 
# 67781   447 
table(dups$case_status,dups$in_prev_analysis,useNA = "ifany")
#       0     1
# 0 66929     2
# 1   852   445
all(dat$in_prev_analysis[idx]==dups$in_prev_analysis) #T
table(dat$j_blad_cancer)
# 0      1 
# 152367   2520 
table(dat$j_blad_cancer_first)
# 0    1 
# 483 2037
table(dat$exc_blad_is_first_dx,dat$j_blad_cancer_first)
table(dat$exc_blad_is_first_dx)
# 0    1 
# 483 2037
table(dat$exc_in_tgwas_population)
# 0      1 
# 44325 110562 
table(dat$in_prev_analysis,dat$exc_blad_is_first_dx,useNA="ifany")
#       0      1    NaN
# 0    473   1551 152246
# 1     10    486    121

# selsamples=dat$plco_id[which(dat$in_prev_analysis==0)]
# selsamples=dups$plco_id[dups$in_prev_analysis==0]
selsamples=dat1$plco_id
selsamples1=paste0(selsamples,"_",selsamples)
write.table(selsamples1,"../result/newPLCA_plinksamples.txt",row.names = F,quote=F)
#check samples in vcffiles
check_vcfsamples=function(infolder="/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/Oncoarray/European/")
{
  vcfsamples=read.table(paste0(infolder,"samples.txt"))
  print(sum(vcfsamples$V1 %in% selsamples1))
  res=intersect(vcfsamples$V1,selsamples1)
  return(res)
}
selsamples1=paste0(dat$plco_id,"_",dat$plco_id)
selsamples1=paste0(dat1$plco_id,"_",dat1$plco_id)
s1=check_vcfsamples()
s2=check_vcfsamples("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/OmniX/European/")
s3=check_vcfsamples("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/Omni25/European/")
s4=check_vcfsamples("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/GSA/batch1/European/")
s5=check_vcfsamples("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/GSA/batch2/European/")
s6=check_vcfsamples("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/GSA/batch3/European/")
s7=check_vcfsamples("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/GSA/batch4/European/")
sum(length(c(s1,s2,s3,s4,s5,s6,s7)))
sum(length(unique(c(s1,s2,s3,s4,s5,s6,s7))))
allvcfsamples=c(s1,s2,s3,s4,s5,s6,s7)
tmp=unlist(strsplit(allvcfsamples,"_"))
allvcfsamples=tmp[seq(1,length(tmp),2)]
idx=match(allvcfsamples,dat$plco_id)
table(dat$in_prev_analysis[idx])
table(dat$exc_blad_is_first_dx[idx],useNA = "ifany")
# 1   NaN 
# 522 40739 
idx=match(allvcfsamples,dups$plco_id)
table(dups$case_status[idx])
idx=match(selsamples,dat$plco_id)
tmp=data.frame(ID=selsamples1)
write.table(tmp,file="../result/newPLCA_bcftoolssamples.txt",row.names = F,col.names=F,quote=F)
table(dat$exc_in_tgwas_population,dat$exc_blad_is_first_dx,useNA="ifany")
dat1=read.csv("../data/package-plco-957/Bladder File/plco_genotype_042022.csv")
table(dat1$plco_id %in% selsamples)
table(dat1$population)
table(dat1$plco_id %in% dat$plco_id) 
# TRUE 
# 41261 
idx=match(dat1$plco_id,dat$plco_id)
table(dat$in_prev_analysis[idx])
table(dat$exc_blad_is_first_dx[idx])
table(dat1$age==dat$age[idx])
table(dat1$exc_blad_is_first_dx,dat$exc_blad_is_first_dx[idx],useNA="ifany")
pheno=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
pheno$TGSID[which(pheno$TGSID=="#N/A")]=NA
table(dat$plco_id %in% pheno$study_pid)

tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/imputation/newPLCO/European/chr22_filtered.psam"))
tmp1=unlist(strsplit(tmp$`#IID`,"_"))
tmp$ID=tmp1[seq(1,length(tmp1),2)]
table(tmp$ID %in% dat$plco_id)

#!/usr/bin/env Rscript
#run gwas on 610K samples (exclude 3 studies)
plink="/usr/local/apps/plink/1.9/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
.libPaths(c("/data/wangx53",.libPaths()))

#change TGSID to study_pid
pheno=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
pheno$TGSID[which(pheno$TGSID=="#N/A")]=NA
tmp=fread("../result/imputation/merged.psam")
changesampleid=function(id=tmp$`#IID`)
{
  dat=data.frame(id=id,study_pid=NA,TGSID=NA)
  tmp=intersect(dat$id,pheno$study_pid)
  idx1=match(tmp,dat$id)
  idx2=match(tmp,pheno$study_pid)
  dat$study_pid[idx1]=pheno$study_pid[idx2]
  tmp=intersect(dat$id,pheno$TGSID)
  idx1=match(tmp,dat$id)
  idx2=match(tmp,pheno$TGSID)
  dat$study_pid[idx1]=pheno$study_pid[idx2]
  if (sum(is.na(dat$study_pid))>0) stop("samples ID missing")
  return(dat$study_pid)
}

library(data.table)
#work on 610K individual data
six10dat=read.table("/data/BB_Bioinformatics/ProjectData/Bladder/Phenotype_data/Overall/610K/snptest.def",header=T)
six10dat=six10dat[-1,]
six10dat$ID_1=changesampleid(id=six10dat$ID_1)
six10dat$ID_2=six10dat$ID_1
#tmp=fread("../result/imputation/chr22.traw",nrows=10000)

#system("cp ../result/imputation/merged.psam ../result/imputation/merged.psam.orig")
#system("cp ../result/imputation/merged.fam ../result/imputation/merged.fam.orig")
tmp=fread("../result/imputation/merged.psam")
tmp$`#IID`=changesampleid(id=tmp$`#IID`)
idx=match(tmp$`#IID`,pheno$study_pid)
tmp$SEX=pheno$gender[idx]
tmp$SEX[which(tmp$SEX=="FEMALE")]=2
tmp$SEX[which(tmp$SEX=="MALE")]=1
tmp$Pheno=1
tmp$Pheno[which(pheno$casecontrol[idx]=="CASE")]=2
write.table(tmp,file="../result/imputation/merged.psam",sep=" ",quote=F,row.names = F)
tmp=fread("../result/imputation/merged.fam")
tmp$V2=changesampleid(id=tmp$V2)
idx=match(tmp$V2,pheno$study_pid)
tmp$V5=pheno$gender[idx]
tmp$V5[which(tmp$SEX=="FEMALE")]=2
tmp$V5[which(tmp$SEX=="MALE")]=1
tmp$V6=1
tmp$V6[which(pheno$casecontrol[idx]=="CASE")]=2
write.table(tmp,file="../result/imputation/merged.fam",sep=" ",quote=F,row.names = F)

#pick samples from 3 studies
idx=which(pheno$study %in% c("CPSII") & pheno$Illumina_Array=="610K")
table(pheno$casecontrol[idx])
CPSIIsample=pheno$study_pid[idx]
idx=which(pheno$study %in% c("NEBL"))
table(pheno$casecontrol[idx])
idx=which(pheno$study %in% c("MDACC") & pheno$Illumina_Array=="610K")
table(pheno$casecontrol[idx])
idx=which(pheno$study %in% c("CPSII","NEBL","MDACC") & pheno$Illumina_Array=="610K")
table(pheno$casecontrol[idx])
all3studies=pheno$study_pid[idx]
all3studies=intersect(all3studies,six10dat$ID_1)
set.seed(1000)
idx=sample(1:length(all3studies),0.5*length(all3studies))
tunsamples=all3studies[idx]
valsamples=all3studies[!all3studies %in% tunsamples]
pheno$cig_cat=factor(pheno$cig_cat,levels=c("NEVER","FORMER","OCCASIONAL","CURRENT"))
pheno$y=0
pheno$y[which(pheno$casecontrol=="CASE")]=1
model1 <- glm(y~cig_cat, data=pheno[pheno$study_pid %in% valsamples,],family = "binomial")
summary(model1)
coeff1 =coef(model1)
exp(coeff1)
exp(confint(model1))

gwassamples=six10dat$ID_1[!six10dat$ID_1 %in% all3studies]
idx=match(gwassamples,pheno$study_pid)
table(pheno$y[idx])
#2650 0, 1024 1
tmp=data.frame(IID=gwassamples)
write.table(tmp,file="../result/GWASsamples.txt",sep=" ",row.names = F,quote=F)
#covariates
idx=match(gwassamples,six10dat$ID_1)
covdat=cbind(tmp,six10dat[idx,4:6])
idx1=match(gwassamples,pheno$study_pid)
table(six10dat$case_control[idx],pheno$y[idx1])
write.table(covdat,file="../result/GWAS.covar",row.names = F,col.names=T,sep=" ",quote=F)
# phenodat=cbind(tmp,case=NA)
# phenodat$case=as.numeric(six10dat$case_control[idx])
# phenodat$case=phenodat$case+1
# write.table(phenodat,file="../result/GWAS.pheno",row.names = F,col.names=T,sep=" ",quote=F)
cmd=paste0(plink2," --pfile ../result/imputation/merged "," --keep ../result/GWASsamples.txt --maf 0.01 --covar  ../result/GWAS.covar --logistic hide-covar --ci 0.95 --out ../result/six10kGWAS")
system(cmd)
tmp1=fread("../result/six10kGWAS.Pheno.glm.logistic.hybrid")
colnames(tmp1)[1:3]=c("CHR","BP","SNP")

tmp1$A2=tmp1$REF
idx=which(tmp1$A1==tmp1$A2)
tmp1$A2[idx]=tmp1$ALT[idx]
tmp1$CASES_Num=1024
tmp1$CONTROLS_Num=2650
colnames(tmp1)[which(colnames(tmp1)=="LOG(OR)_SE")]="SE"
cmd=paste0(plink2," --pfile ../result/imputation/merged "," --keep ../result/GWASsamples.txt --freq --out ../result/six10kGWAS")
system(cmd)
freqdat=fread("../result/six10kGWAS.afreq")
idx=match(tmp1$SNP,freqdat$ID)
freqdat=freqdat[idx,]
tmp1$A1_EAF=freqdat$ALT_FREQS
idx=which(freqdat$ALT!=tmp1$A1)
tmp1$A1_EAF[idx]=1-tmp1$A1_EAF[idx]
fwrite(tmp1,file="../result/six10kGWAS.txt",row.names = F,sep="\t",quote=F)
onem=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/1M/out.plink.rsq03.01.txt"))

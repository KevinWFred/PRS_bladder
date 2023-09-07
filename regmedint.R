#!/usr/bin/env Rscript
#mediator smoking on PRS

library(data.table)
library(regmedint)
library(dplyr)
#tuning samples
famtrain=read.table("../result/six10k_train.fam")
#validation samples
famtest=read.table("../result/six10k_test.fam")
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

phenotype=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
pheno=phenotype
idx=match(c(famtrain$V2,famtest$V2),phenotype$TGSID)
idx1=which(is.na(idx))
idx2=match(c(famtrain$V2,famtest$V2)[idx1],phenotype$study_pid)
idx[idx1]=idx2
sum(is.na(idx))
phenotype=phenotype[idx,]
phenotype$ID=c(famtrain$V2,famtest$V2) #ID for all samples
phenotype$y=0
phenotype$y[which(phenotype$casecontrol=="CASE")]=1
pcadat=fread("../result/merged.eigenvec")
colnames(pcadat)[2:ncol(pcadat)]=c("ID",paste0("EV",1:20))
tmp=changesampleid(id=pcadat$ID)
idx=match(phenotype$ID,tmp)
phenotype=cbind(phenotype,pcadat[idx,3:12])
idx=which(phenotype=="MISSING",arr.ind=T)
phenotype[idx]=NA
table(phenotype$cig_cat)
# CURRENT     FORMER      NEVER OCCASIONAL 
# 874       2695       1669         84
phenotype$type=NA
phenotype$type[which(phenotype$ID %in% famtrain$V2)]="tuning"
phenotype$type[which(phenotype$ID %in% famtest$V2)]="validation"
phenotype=phenotype[phenotype$cig_cat %in% c("NEVER","FORMER","CURRENT"),]
phenotype$cig_cat=factor(phenotype$cig_cat,levels=c("NEVER","FORMER","CURRENT"))
phenotype$Age=as.numeric(phenotype$Age)
phenotype0=phenotype
smoke_bin = model.matrix(~as.factor(cig_cat), data = phenotype)[,-1]
colnames(smoke_bin) = c("Former", "Current")
rownames(smoke_bin)=phenotype$ID
phenotype=phenotype[,c("ID","y","Age","gender","type",paste0("EV",1:10))]
phenotype=cbind(phenotype,smoke_bin)
load("../result/PRS_res.RData")
allprs=GWAS24res$allprs
data_com = inner_join(phenotype,allprs, by = c("ID"))
data_com_control = data_com[data_com$y==0,]
mean_prs = mean(data_com_control$prs,na.rm = T)
se_prs = sd(data_com_control$prs,na.rm = T)
data_com$prs_sd = (data_com$prs-mean_prs)/se_prs
data_com$gender=as.numeric(factor(data_com$gender))-1
mean_EV=rep(NA,10)
for (i in 1:10)
{
  mean_EV[i]=mean(data_com[,paste0("EV",i)])
}
#FORMER
C = c(mean(data_com$Age,na.rm=T), 0,0,mean_EV)

result <- regmedint(data = data_com,
                    ## Variables
                    yvar = "y",
                    avar = "prs_sd",
                    mvar = "Former",
                    cvar = c("Age", "gender", "Current",paste0("EV",1:10)),
                    ## Values at which effects are evaluated
                    a0 = 0,
                    a1 = 1,
                    m_cde = 0,
                    c_cond = C,
                    ## Model types
                    mreg = "logistic",
                    yreg = "logistic",
                    ## Additional specification
                    interaction = F,
                    casecontrol = F)
summary(result)
#In the result
#row 2 pnde is the direct effect
#row 3 tnie is the indirect effect
#row 6 te is the total effect
#row 7 pm is the proportion of indirect effect over total effect
#          est         se         Z            p      lower      upper
# cde  0.416929098 0.029780778 13.999940 0.0000000  0.358559845 0.47529835
# pnde 0.416929098 0.029780778 13.999940 0.0000000  0.358559845 0.47529835
# tnie 0.008444875 0.005302697  1.592562 0.1112584 -0.001948221 0.01883797
# tnde 0.416929098 0.029780778 13.999940 0.0000000  0.358559845 0.47529835
# pnie 0.008444875 0.005302697  1.592562 0.1112584 -0.001948221 0.01883797
# te   0.425373973 0.030277204 14.049315 0.0000000  0.366031743 0.48471620
# pm   0.024271089 0.014988073  1.619360 0.1053698 -0.005104994 0.05364717

#Current
C = c(mean(data_com$Age,na.rm=T), 0,0,mean_EV)

result1 <- regmedint(data = data_com,
                    ## Variables
                    yvar = "y",
                    avar = "prs_sd",
                    mvar = "Current",
                    cvar = c("Age", "gender", "Former",paste0("EV",1:10)),
                    ## Values at which effects are evaluated
                    a0 = 0,
                    a1 = 1,
                    m_cde = 0,
                    c_cond = C,
                    ## Model types
                    mreg = "logistic",
                    yreg = "logistic",
                    ## Additional specification
                    interaction = F,
                    casecontrol = F)
summary(result1)
#          est         se         Z            p      lower      upper
# cde  0.41692910 0.02978078 13.999940 0.000000e+00 0.35855984 0.47529835
# pnde 0.41692910 0.02978078 13.999940 0.000000e+00 0.35855984 0.47529835
# tnie 0.06527649 0.01583737  4.121674 3.761289e-05 0.03423581 0.09631716
# tnde 0.41692910 0.02978078 13.999940 0.000000e+00 0.35855984 0.47529835
# pnie 0.06527649 0.01583737  4.121674 3.761289e-05 0.03423581 0.09631716
# te   0.48220558 0.03385373 14.243793 0.000000e+00 0.41585349 0.54855768
# pm   0.16517227 0.03539015  4.667182 3.053582e-06 0.09580886 0.23453568
#!/usr/bin/env Rscript
#mediator smoking on PRS
.libPaths(c("/data/wangx53",.libPaths()))

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
#allprs=LDpredres$allprs
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
C = c(mean(data_com$Age,na.rm=T),0,mean_EV)

result <- regmedint(data = data_com,
                    ## Variables
                    yvar = "y",
                    avar = "prs_sd",
                    mvar = "Former",
                    cvar = c("Age", "gender",paste0("EV",1:10)),
                    ## Values at which effects are evaluated
                    a0 = 0,
                    a1 = 1,
                    m_cde = 0,
                    c_cond = C,
                    ## Model types
                    mreg = "logistic",
                    yreg = "logistic",
                    ## Additional specification
                    interaction = T,
                    casecontrol = T)
tmp=summary(result)
#In the result
#row 2 pnde is the direct effect
#row 3 tnie is the indirect effect
#row 6 te is the total effect
#row 7 pm is the proportion of indirect effect over total effect
round(tmp$summary_myreg,3)
tmp=summary(result)
write.csv(tmp$summary_myreg,file="../result/regmedint_formersmoker_mediator.csv",row.names = T,quote=F)
#Current
C = c(mean(data_com$Age,na.rm=T), 0,mean_EV)

result1 <- regmedint(data = data_com,
                    ## Variables
                    yvar = "y",
                    avar = "prs_sd",
                    mvar = "Current",
                    cvar = c("Age", "gender", paste0("EV",1:10)),
                    ## Values at which effects are evaluated
                    a0 = 0,
                    a1 = 1,
                    m_cde = 0,
                    c_cond = C,
                    ## Model types
                    mreg = "logistic",
                    yreg = "logistic",
                    ## Additional specification
                    interaction = T,
                    casecontrol = T)
summary(result1)
#          est         se         Z            p      lower      upper
# cde   0.41486325 0.02947408 14.075528 0.0000000  0.35709512 0.472631388
# pnde  0.41486325 0.02947408 14.075528 0.0000000  0.35709512 0.472631388
# tnie -0.01021555 0.00803905 -1.270741 0.2038208 -0.02597180 0.005540697
# tnde  0.41486325 0.02947408 14.075528 0.0000000  0.35709512 0.472631388
# pnie -0.01021555 0.00803905 -1.270741 0.2038208 -0.02597180 0.005540697
# te    0.40464770 0.03052789 13.255018 0.0000000  0.34481414 0.464481261
# pm   -0.03085419 0.02496329 -1.235982 0.2164652 -0.07978135 0.018072971
tmp=summary(result1)
write.csv(tmp$summary_myreg,file="../result/regmedint_currentsmoker_mediator.csv",row.names = T,quote=F)

# #Run smoking-->PRS-->outcome
# mean_EV=rep(NA,10)
# for (i in 1:10)
# {
#   mean_EV[i]=mean(data_com[,paste0("EV",i)])
# }
# #FORMER
# C = c(mean(data_com$Age,na.rm=T),0,mean_EV)
# 
# result <- regmedint(data = data_com,
#                     ## Variables
#                     yvar = "y",
#                     avar = "Former",#"prs_sd",
#                     mvar = "prs_sd",#"Former",
#                     cvar = c("Age", "gender",paste0("EV",1:10)),
#                     ## Values at which effects are evaluated
#                     a0 = 0,
#                     a1 = 1,
#                     m_cde = 0,
#                     c_cond = C,
#                     ## Model types
#                     mreg = "linear",
#                     yreg = "logistic",
#                     ## Additional specification
#                     interaction = T,
#                     casecontrol = F)
# summary(result)
# #In the result
# #row 2 pnde is the direct effect
# #row 3 tnie is the indirect effect
# #row 6 te is the total effect
# #row 7 pm is the proportion of indirect effect over total effect
# # #          est         se         Z            p      lower      upper
# # cde   0.170500310 0.06029951  2.8275570 0.004690467  0.052315435 0.28868519
# # pnde  0.133767168 0.06389811  2.0934449 0.036309458  0.008529178 0.25900516
# # tnie -0.004651169 0.01162840 -0.3999838 0.689168472 -0.027442405 0.01814007
# # tnde  0.134438699 0.06364375  2.1123630 0.034655327  0.009699244 0.25917815
# # pnie -0.005322700 0.01330521 -0.4000464 0.689122368 -0.031400429 0.02075503
# # te    0.129115999 0.06493501  1.9883882 0.046768775  0.001845725 0.25638627
# # pm   -0.038488238 0.10126085 -0.3800900 0.703878588 -0.236955849 0.15997937
# tmp=summary(result)
# write.csv(tmp$summary_myreg,file="../result/regmedint_prs_mediate_formersmoker.csv",row.names = T,quote=F)
# 
# #Current
# C = c(mean(data_com$Age,na.rm=T), 0,mean_EV)
# 
# result1 <- regmedint(data = data_com,
#                      ## Variables
#                      yvar = "y",
#                      avar = "Current",#"prs_sd",
#                      mvar = "prs_sd",#"Current",
#                      cvar = c("Age", "gender", paste0("EV",1:10)),
#                      ## Values at which effects are evaluated
#                      a0 = 0,
#                      a1 = 1,
#                      m_cde = 0,
#                      c_cond = C,
#                      ## Model types
#                      mreg = "linear",
#                      yreg = "logistic",
#                      ## Additional specification
#                      interaction = T,
#                      casecontrol = F)
# summary(result1)
# #         est         se         Z            p      lower      upper
# # cde  1.11887311 0.08653852 12.929192 0.000000e+00 0.94926072 1.28848549
# # pnde 1.26791092 0.10173358 12.463052 0.000000e+00 1.06851676 1.46730508
# # tnie 0.09649685 0.02661908  3.625101 2.888487e-04 0.04432441 0.14866929
# # tnde 1.30283021 0.10984481 11.860645 0.000000e+00 1.08753835 1.51812208
# # pnie 0.06157756 0.01573569  3.913241 9.106548e-05 0.03073617 0.09241895
# # te   1.36440777 0.11120746 12.269031 0.000000e+00 1.14644516 1.58237039
# # pm   0.12356103 0.03058311  4.040172 5.341208e-05 0.06361923 0.18350284
# 
# tmp=summary(result1)
# write.csv(tmp$summary_myreg,file="../result/regmedint_prs_mediate_currentsmoker.csv",row.names = T,quote=F)

#each of 24 snps to be mediator

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
#phenotype=phenotype[,c("ID","y","Age","gender","type",paste0("EV",1:10))]
phenotype=cbind(phenotype,smoke_bin)

gwasdat=read.table("../result/gwas_24snp/Bladder_gwas24.traw",header=T)
rownames(gwasdat)=gwasdat$SNP
gwasdat=gwasdat[,7:ncol(gwasdat)]
colnames(gwasdat)=gsub("X0_","",colnames(gwasdat))
table(phenotype$study_pid %in% colnames(gwasdat))
idx=which(!phenotype$study_pid %in% colnames(gwasdat))
table(phenotype$TGSID %in% colnames(gwasdat))
#change gwasdat colnames to study_pid
tmp=intersect(phenotype$TGSID, colnames(gwasdat))
idx1=match(tmp,phenotype$TGSID)
idx2=match(tmp,colnames(gwasdat))
colnames(gwasdat)[idx2]=phenotype$study_pid[idx1]
table(phenotype$study_pid %in% colnames(gwasdat))
tmp=intersect(phenotype$study_pid,colnames(gwasdat))
idx1=match(tmp,phenotype$study_pid)
idx2=match(tmp,colnames(gwasdat))
gwasdat=gwasdat[,idx2]
phenotype=phenotype[idx1,]
all(phenotype$ID == colnames(gwasdat))

data=phenotype
data$gender=as.numeric(factor(data$gender))-1
mean_EV=rep(NA,10)
for (i in 1:10)
{
  mean_EV[i]=mean(data[,paste0("EV",i)])
}

C = c(mean(data$Age,na.rm=T),0,mean_EV)

run_pregmedint=function(i=1)
{
  idx=match(colnames(gwasdat),data$study_pid)
  data_com=cbind.data.frame(data,snp=unlist(gwasdat[i,idx]))
  data_com_control = data_com[data_com$y==0,]
  mean_snp = mean(data_com_control$snp,na.rm = T)
  se_snp = sd(data_com_control$snp,na.rm = T)
  data_com$snp_sd = (data_com$snp-mean_snp)/se_snp
  result <- regmedint(data = data_com,
                      ## Variables
                      yvar = "y",
                      avar = "snp_sd",
                      mvar = "Former",
                      cvar = c("Age", "gender",paste0("EV",1:10)),
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = C,
                      ## Model types
                      mreg = "logistic",
                      yreg = "logistic",
                      ## Additional specification
                      interaction = T,
                      casecontrol = T)
  tmp=summary(result)
  resformer=tmp$summary_myreg
  
  result1 <- regmedint(data = data_com,
                       ## Variables
                       yvar = "y",
                       avar = "snp_sd",
                       mvar = "Current",
                       cvar = c("Age", "gender", paste0("EV",1:10)),
                       ## Values at which effects are evaluated
                       a0 = 0,
                       a1 = 1,
                       m_cde = 0,
                       c_cond = C,
                       ## Model types
                       mreg = "logistic",
                       yreg = "logistic",
                       ## Additional specification
                       interaction = T,
                       casecontrol = T)
  tmp=summary(result1)
  rescurrent=tmp$summary_myreg
  return(list(resformer=resformer,rescurrent=rescurrent))
}
resformer=rescurrent=list()
for (i in 1:24)
{
  res=run_pregmedint(i=i)
  resformer[[rownames(gwasdat)[i]]]=res$resformer
  rescurrent[[rownames(gwasdat)[i]]]=res$rescurrent
}
allres=data.frame(snp=rownames(gwasdat),Former_pnde_est=NA,Formerpnde_p=NA,
                  Former_tnie_est=NA,Former_tnie_p=NA,
                  Former_te_est=NA,Former_te_p=NA,
                  Former_pm_est=NA,Former_pm_p=NA,
                  Current_pnde_est=NA,Currentpnde_p=NA,
                  Current_tnie_est=NA,Current_tnie_p=NA,
                  Current_te_est=NA,Current_te_p=NA,
                  Current_pm_est=NA,Current_pm_p=NA)
for(i in 1:24)
{
  tmp=resformer[[rownames(gwasdat)[i]]]
  tmp1=rescurrent[[rownames(gwasdat)[i]]]
  allres[i,2:9]=c(tmp[2,c(1,4)],tmp[3,c(1,4)],tmp[6,c(1,4)],tmp[7,c(1,4)])
  allres[i,10:ncol(allres)]=c(tmp1[2,c(1,4)],tmp1[3,c(1,4)],tmp1[6,c(1,4)],tmp1[7,c(1,4)])
}
write.csv(allres[order(allres$Current_tnie_p),],file="../result/regmedint_24gwas.csv",row.names = F,quote=F)

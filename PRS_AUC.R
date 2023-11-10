#!/usr/bin/env Rscript
#to get all PRS and AUC for FLCCA samples
#Adjsted AUC and AUC for PRS
.libPaths(c("/data/wangx53",.libPaths()))
library("readxl")
library(data.table)
library(pROC)
#library(SuperLearner)
library(RISCA)
library(boot)
library(dbplyr)

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
phenotype$cig_cat=factor(phenotype$cig_cat,levels=c("NEVER","FORMER","OCCASIONAL","CURRENT"))


RISCA_AUC=function(allprs=CT_never$allprs)
{
  phenotype2=phenotype[match(famtest$V2,phenotype$ID),]
  pheno.prs=merge(phenotype2,allprs,by="ID")
  pheno.prs=myscale(pheno.prs)
  
  # roc_obj_pv = roc.binary(status = "y", estimator = "pv",
  #                         variable = "prs",
  #                         confounders = ~EV1+EV2+EV3+EV4+EV5+EV6+
  #                           EV7+EV8+EV9+EV10+AGE,
  #                         data = pheno.prs,
  #                         precision=seq(0.05,0.95, by=0.05))
  confounders = c(paste0("EV",1:10),"Age_cat","gender")
  idx=match(confounders,colnames(pheno.prs))
  idx=complete.cases(pheno.prs[,idx])
  pheno.prs=pheno.prs[idx,]
  roc_obj_ipw = roc.binary(status = "y", estimator = "ipw",
                           variable = "prs",
                           #confounders = ~1,
                           confounders = ~EV1+EV2+EV3+EV4+EV5+EV6+
                             EV7+EV8+EV9+EV10+Age_cat+gender,
                           data = pheno.prs,
                           precision=seq(0.05,0.95, by=0.05))
  #print(roc_obj_pv$auc)
  #print(roc_obj_ipw$auc)
  return(roc_obj_ipw$auc)
}

AUCBoot = function(data,indices){
  boot_data = data[indices, ]
  model4 <- glm(y~prs, data=boot_data,family = "binomial")
  predicted4 <- predict(model4,boot_data, type="response")
  auc=as.numeric(pROC::auc(boot_data$y,predicted4,quiet=T))
  return(c(auc))
}

AUCadjBoot = function(data,indices){
  boot_data = data[indices, ]
  confounders = c(paste0("EV",1:10),"Age_cat","gender")
  idx=match(confounders,colnames(data))
  idx=complete.cases(boot_data[,idx])
  boot_data=boot_data[idx,]
  roc_obj_ipw = roc.binary(status = "y", estimator = "ipw",
                           variable = "prs",
                           confounders = ~EV1+EV2+EV3+EV4+EV5+EV6+
                             EV7+EV8+EV9+EV10+Age_cat+gender,
                           data = boot_data,
                           precision=seq(0.05,0.95, by=0.05))
  auc=roc_obj_ipw$auc
  return(c(auc))
}

#standardize based on controls
myscale=function(pheno.prs)
{
  idx=pheno.prs$casecontrol=="CONTROL"
  pheno.prs$prs=(pheno.prs$prs-mean(pheno.prs$prs[idx]))/sd(pheno.prs$prs[idx])
  return(pheno.prs)
}

get_aucs=function(allprs)
{
  auc=data.frame(tun=0,val=0)
  pheno.prs=merge(phenotype,allprs,by="ID")
  model1 <- glm(y~prs, data=pheno.prs[pheno.prs$ID %in% famtrain$V2,],family = "binomial")
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtrain$V2,], type="response")
  auc$tun=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtrain$V2],predicted1,quiet=T))
  model1 <- glm(y~prs, data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$val=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  
  print("auc:")
  print(auc)
  auc_val=data.frame(auc_low=0,auc_high=0)
  #on validation
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc_val$auc_low=tmp$bca[4]
  auc_val$auc_high=tmp$bca[5]
  print(auc_val)
  
  #only for validation
  aucadj=data.frame(auc=0,auc_low=0,auc_high=0)
  aucadj$auc=RISCA_AUC(allprs=allprs)
  aucadj_val=data.frame(auc_low=0,auc_high=0)
  boot_aucadj = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCadjBoot, R = 10000)
  tmp=boot.ci(boot_aucadj,type="bca")
  aucadj$auc_low=tmp$bca[4]
  aucadj$auc_high=tmp$bca[5]
  print(aucadj)
  
  return(list(auc=auc,auc_val=auc_val,aucadj=aucadj))
}

prs_or_quantiles=function(allprs=CTprs$allprs,main="",ytext=0.25)
{
  phenotype2=phenotype[match(famtest$V2,phenotype$ID),]
  pheno.prs=merge(phenotype2,allprs,by="ID")
  
  #pheno.prs$prs=scale(pheno.prs$prs)
  pheno.prs=myscale(pheno.prs)
  model1 <- glm(y~prs+Age_cat+gender+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno.prs,family = "binomial")
  coeff1 =coef(model1)
  OR=data.frame(or=0,or_low=0,or_high=0)
  OR$or=exp(coeff1)[2]
  tmp1=confint(model1)
  OR$or_low=exp(tmp1[2,1])
  OR$or_high=exp(tmp1[2,2])
  print(OR)
  idx=pheno.prs$casecontrol=="CONTROL"
  quantiles=quantile(pheno.prs$prs[idx],c(0,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,1))
  pheno.prs$prsquantile=cut(pheno.prs$prs,quantiles,labels=c("< 5%","5-10%","10-20%","20-40%","40-60%","60-80%","80-90%","90-95%","> 95%"),include.lowest = T)
  pheno.prs$prsquantile=factor(pheno.prs$prsquantile,levels=c("40-60%","< 5%","5-10%","10-20%","20-40%","60-80%","80-90%","90-95%","> 95%"))
  idx1=which(pheno.prs$prs>max(pheno.prs$prs[idx]))
  if (length(idx1)>0)
  {
    pheno.prs$prsquantile[idx1]=factor("> 95%",levels=c("40-60%","< 5%","5-10%","10-20%","20-40%","60-80%","80-90%","90-95%","> 95%"))
  }
  idx1=which(pheno.prs$prs<min(pheno.prs$prs[idx]))
  if (length(idx1)>0)
  {
    pheno.prs$prsquantile[idx1]=factor("< 5%",levels=c("40-60%","< 5%","5-10%","10-20%","20-40%","60-80%","80-90%","90-95%","> 95%"))
  }
  model2 <- glm(y~prsquantile+Age_cat+gender+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno.prs,family = "binomial")
  tmp=summary(model2)$coefficient
  tmp1=confint(model2)
  n=length(unique(pheno.prs$prsquantile))
  n1=which(rownames(tmp1)=="prsquantile20-40%")
  if (!any(is.na(tmp1[2:n,])))
  {
    OR_quantile=data.frame(or=c(exp(tmp[2:n1,1]),1,exp(tmp[(n1+1):n,1])),or_low=c(exp(tmp1[2:n1,1]),1,exp(tmp1[(n1+1):n,1])),
                           or_high=c(exp(tmp1[2:n1,2]),1,exp(tmp1[(n1+1):n,2])))
    ymin=max(min(OR_quantile$or_low) -(max(OR_quantile$or_high)-min(OR_quantile$or_low))*0.02,0.1)
    ymax=min(max(OR_quantile$or_high) +(max(OR_quantile$or_high)-min(OR_quantile$or_low))*0.1,30)
    par(mar = c(4.5, 4.5, 2, 2) + 0.2)
    plot(1:n,OR_quantile$or,xlab="",ylab="Odds ratio (95% CI)",ylim=c(ymin,ymax),cex.lab=1.2,cex.axis=1.2,pch=20,col="blue",xaxt="n",yaxt="n",log="y",main=main)
    for (i in 1:n)
    {
      lines(c(i,i),c(OR_quantile$or_low[i],OR_quantile$or_high[i]),lwd=2,col="blue")
      if (OR_quantile$or_low[i]<OR_quantile$or_high[i])
      {
        lines(c(i-0.1,i+0.1),c(OR_quantile$or_low[i],OR_quantile$or_low[i]),lwd=2,col="blue")
        lines(c(i-0.1,i+0.1),c(OR_quantile$or_high[i],OR_quantile$or_high[i]),lwd=2,col="blue")
      }
    }
  }
  axis(side=1,at=1:n,labels=F)
  
  #axis(side=1,at=1:n,c("< 1%","1-5%","5-10%","10-20%","20-40%","40-60%","60-80%","80-90%","90-95%","95-99%","> 99%"),las=2)
  axis(side=2,at=c(seq(0,1,0.2),seq(2,20,2)),las=2)
  for ( i in c(seq(0.2,1,0.2),2,4,8,16))
    abline(h=i,lty=3,col="gray")
  # text(x=1:n, y=rep(-0.5,n),
  #      labels=c("< 1%","1-5%","5-10%","10-20%","20-40%","40-60%","60-80%","80-90%","90-95%","95-99%","> 99%"),  xpd=NA)
  mtext("PRS percentile",side = 1,line=3.5,cex=1.2)
  text(x = 1:n,
       ## Move labels to just below bottom of chart.
       y = ytext,
       ## Use names from the data list.
       labels = c("< 5%","5-10%","10-20%","20-40%","40-60%","60-80%","80-90%","90-95%","> 95%"),
       ## Change the clipping region.
       xpd =T,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       
       ## Adjust the labels to almost 100% right-justified.
       adj = 1,
       ## Increase label size.
       cex = 1)
  return(OR)
}

#CT
#p-value cutoffs (consistent with ../result/range_list)
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 
CT_prs=function(prsprefix_tun="../result/six10k_tun.p_value_",prsprefix_val="../result/six10k_val.p_value_")
{
  #find the optimal cutoff
  auc_tun=rep(0,length(pthres))
  for (i in 1:length(pthres))
  {
    prs=read.table(paste0(prsprefix_tun,i,".sscore"))
    colnames(prs)[2]="ID"
    colnames(prs)[6]="prs"
    pheno.prs=merge(phenotype,prs,by="ID")
    
    model1 <- glm(y~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    auc_tun[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
   
  }
  idx_optimal=which.max(auc_tun)
  clumpedsnp=fread("../result/six10k.clumped")
  sum(clumpedsnp$P<=pthres[idx_optimal]) #20 snps (P<5e-7) used
  CTselsnps=clumpedsnp$SNP[which(clumpedsnp$P<=pthres[idx_optimal])]
  write.table(CTselsnps,file="../result/CTselsnps.txt",row.names = F,col.names = F,quote=F)
  prs_tun=read.table(paste0(prsprefix_tun,idx_optimal,".sscore"))
  colnames(prs_tun)[2]="ID"
  colnames(prs_tun)[6]="prs"
  prs_tun=prs_tun[,c(2,6)]
  
  prs=read.table(paste0(prsprefix_val,idx_optimal,".sscore"))
  colnames(prs)[2]="ID"
  colnames(prs)[6]="prs"
  prs_val=prs[,c(2,6)]
  allprs=rbind(prs_tun,prs_val)
  
  aucs=get_aucs(allprs=allprs)
  
  return(list(aucs=aucs,idx_optimal=idx_optimal,allprs=allprs))
}
CTres=CT_prs()
# $auc
# tun       val
# 1 0.5831845 0.5991669
# 
# $auc_val
# auc_low  auc_high
# 1 0.5778936 0.6197845
# 
# $aucadj
# auc   auc_low  auc_high
# 1 0.6007601 0.5795602 0.6225562
png("../result/CT_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
CTres_OR=prs_or_quantiles(allprs=CTres$allprs,main="CT",ytext = 0.12)
dev.off()
# or   or_low  or_high
# 1  1.451246 1.340657 1.572681

CT_5e8prs=function(prsprefix_tun="../result/six10k_tun.p_value_",prsprefix_val="../result/six10k_val.p_value_")
{
  
  idx_optimal=1 #5e-8
  clumpedsnp=fread("../result/six10k.clumped")
  sum(clumpedsnp$P<=pthres[1]) #14 snps (P<5e-8) used
  prs_tun=read.table(paste0(prsprefix_tun,idx_optimal,".sscore"))
  colnames(prs_tun)[2]="ID"
  colnames(prs_tun)[6]="prs"
  prs_tun=prs_tun[,c(2,6)]
  
  prs=read.table(paste0(prsprefix_val,idx_optimal,".sscore"))
  colnames(prs)[2]="ID"
  colnames(prs)[6]="prs"
  prs_val=prs[,c(2,6)]
  allprs=rbind(prs_tun,prs_val)
  
  aucs=get_aucs(allprs=allprs)
  
  return(list(aucs=aucs,idx_optimal=idx_optimal,allprs=allprs))
}
CT5e8res=CT_5e8prs()
# $auc
# tun       val
# 1 0.5804069 0.5871434
# 
# $auc_val
# auc_low  auc_high
# 1 0.5649726 0.6082735
# 
# $aucadj
# auc   auc_low  auc_high
# 1 0.5882159 0.5669109 0.6096396

CT5e8res_OR=prs_or_quantiles(allprs=CT5e8res$allprs,main="CT(5e-8)",ytext = 0.27)

#LDpred2----------------------------
LDpred_prs=function(LDpredrda="../result/LDpred_six10k_pred.RData")
{
  load(LDpredrda)
  dim(df_beta)# 997835 snps
  auc_tun=rep(0,ncol(pred_grid))
  phenotype1=phenotype[match(famtrain$V2,phenotype$ID),]
  for (i in 1:ncol(pred_grid))
  {
    pheno.prs=cbind.data.frame(phenotype1,prs=pred_grid[,i])
    model1 <- glm(y~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    auc_tun[i]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  }
  
  idx_optimal=which.max(auc_tun)
  prs_tun=data.frame(ID=names(pred_grid[,idx_optimal]),prs=pred_grid[,idx_optimal])
  prs_val=data.frame(ID=names(pred_grid_val[,idx_optimal]),prs=pred_grid_val[,idx_optimal])
  allprs=rbind(prs_tun,prs_val)
  param_optimal=grid.param[idx_optimal,]
  aucs=get_aucs(allprs=allprs)
  
  return(list(auc=aucs,idx_optimal=idx_optimal,param_optimal=param_optimal,
              allprs=allprs))
}

LDpredres=LDpred_prs()
# $auc
# tun       val
# 1 0.5964123 0.6006132
# 
# $auc_val
# auc_low auc_high
# 1 0.5796134 0.621945
# 
# $aucadj
# auc   auc_low  auc_high
# 1 0.6018158 0.5797877 0.6229344

#write.table(LDpredres$allprs,file="../result/LDpred_prs.csv",sep=",",quote=F,row.names = F)
png("../result/LDpred_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
LDpredres_OR=prs_or_quantiles(allprs=LDpredres$allprs,main="LDpred2",ytext=0.22)
dev.off()

#Lassosum2
Lassosum_prs=function(Lassosumrda="../result/Lassosum2_six10k_pred.RData")
{
  load(Lassosumrda)
  auc_tun=rep(0,ncol(pred_grid))
  phenotype1=phenotype[match(famtrain$V2,phenotype$ID),]
  for (i in 1:ncol(pred_grid))
  {
    if(all(is.na(pred_grid[,i]))) next
    pheno.prs=cbind.data.frame(phenotype1,prs=pred_grid[,i])
    model1 <- glm(y~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    auc_tun[i]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  }
  
  idx_optimal=which.max(auc_tun)
  prs_tun=data.frame(ID=names(pred_grid[,idx_optimal]),prs=pred_grid[,idx_optimal])
  prs_val=data.frame(ID=names(pred_grid_val[,idx_optimal]),prs=pred_grid_val[,idx_optimal])
  allprs=rbind(prs_tun,prs_val)
  param_optimal=grid.param[idx_optimal,]
  aucs=get_aucs(allprs=allprs)
  
  return(list(auc=aucs,idx_optimal=idx_optimal,param_optimal=param_optimal,
              allprs=allprs))
}

Lassosumres=Lassosum_prs()
# $auc
# tun      val
# 1 0.5956919 0.596923
# 
# $auc_val
# auc_low  auc_high
# 1 0.5759035 0.6180192
# 
# $aucadj
# auc   auc_low  auc_high
# 1 0.598126 0.5773293 0.6195446

#write.table(Lassosumres$allprs,file="../result/Lassosum_prs.csv",sep=",",quote=F,row.names = F)
png("../result/Lassosum_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
Lassosumres_OR=prs_or_quantiles(allprs=Lassosumres$allprs,main="Lassosum2",ytext=0.16)
dev.off()
# or or_low or_high
# 1 1.395983 1.3114 1.48701

#PRS-CS
phi = c("phi1e+00","phi1e-02","phi1e-04","phi1e-06")
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

get_prscs=function(prefix="../result/prscs/six10k_pst_eff_a1_b0.5_",bimprefix="../result/six10k_train",outprefix="../result/six10k_train_prscs")
{
  #if (!file.exists(paste0(outprefix,"PRS.sscore")))
  {
    #read posterior SNP effect size
    betadat=NULL
    for (chr in 1:22)
    {
      datlist=list()
      for (i in 1:4)
      {
        myphi=phi[i]
        myfile=paste0(prefix,myphi,"_chr",chr,".txt")
        if (file.exists(myfile))
        {
          datlist[[i]]=fread(myfile)
          
        }else
        {
          warning(paste0(myfile," is not existed"))
        }
      }
      #all(datlist[[1]]$V4==datlist[[4]]$V4)
      betadat_chr=cbind.data.frame(SNP=datlist[[1]]$V2,A1=datlist[[1]]$V4,phi1=datlist[[1]]$V6,phi2=datlist[[2]]$V6,phi4=datlist[[3]]$V6,phi6=datlist[[4]]$V6)
      betadat=rbind(betadat,betadat_chr)
    }
    betaoutfile=paste0(prefix,".prs_coeff")
    write.table(betadat,file=betaoutfile,row.names = F,sep="\t",quote=F)
    
    #use plink to get PRS
    cmd=paste0(plink2," --pfile ",bimprefix," --score-col-nums 3-",ncol(betadat),
               " --score ",betaoutfile," cols=+scoresums,-scoreavgs header no-mean-imputation ",
               " --out ",outprefix,"PRS")
    
    system(cmd)
  }
  
  prsdat=read.table(paste0(outprefix,"PRS.sscore"))
  prsdat=prsdat %>% dplyr::select(V1,V5,V6,V7,V8) %>%
    rename(sample=V1,phi1=V5,phi2=V6,phi4=V7,phi6=V8)
  return(prsdat)
}

PRScs_tun=get_prscs()
PRScs_val=get_prscs(bimprefix="../result/six10k_test",outprefix="../result/six10k_test_prscs")

PRSCS_prs=function(PRS_tun=PRScs_tun,PRS_val=PRScs_val)
{
  auc_tun=data.frame(EUR=rep(0,4))
  auc_val=data.frame(EUR=0,EUR_low=0,EUR_high=0)
  
  phenotype1=phenotype[match(PRS_tun$sample,phenotype$ID),]
  for(i in 1:4)
  {
    pheno.prs=cbind.data.frame(phenotype1,prseur=PRS_tun[,i+1])
    model1 <- glm(y~prseur, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    auc_tun$EUR[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  }
  allbeta=as.data.frame(fread("../result/prscs/six10k_pst_eff_a1_b0.5_.prs_coeff"))
  dim(allbeta) #1053930
  idx_optimal=which.max(auc_tun$EUR)
  allprs=data.frame(ID=c(PRS_tun$sample,PRS_val$sample),prs=c(PRS_tun[,idx_optimal+1],
               PRS_val[,idx_optimal+1]))
  aucs=get_aucs(allprs=allprs)
  return(list(aucs=aucs,allprs=allprs))
}
PRScsres=PRSCS_prs()
# $auc
# tun       val
# 1 0.5809478 0.5981579
# 
# $auc_val
# auc_low  auc_high
# 1 0.5768777 0.6186473
# 
# $aucadj
# auc   auc_low  auc_high
# 1 0.5998702 0.5789898 0.6212744
png("../result/PRScs_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
PRScs_OR=prs_or_quantiles(allprs=PRScsres$allprs,main="PRS-CS",ytext=0.22)
dev.off()
# or   or_low  or_high
# 1 1.334052 1.254841 1.419003
write.table(PRScsres$allprs,file="../result/PRScs_prs.csv",sep=",",quote=F,row.names = F)


#work on 24 GWAS snps
gwas24=function()
{
  prs=read.table("../result/Bladder_gwas24.sscore")
  colnames(prs)[1]="ID"
  colnames(prs)[ncol(prs)]="prs"
  prs$ID=changesampleid(id=prs$ID)
  pheno.prs=merge(phenotype,prs,by="ID")
  allprs=data.frame(ID=pheno.prs$ID,prs=pheno.prs$prs)
  aucs=get_aucs(allprs=allprs)
  
  return(list(aucs=aucs,allprs=allprs))
}
GWAS24res=gwas24()

png("../result/GWAS24_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
GWAS24_OR=prs_or_quantiles(allprs=GWAS24res$allprs,main="GWAS24",ytext = 0.07)
dev.off()

#add smoking into the model
AUCSmokingBoot = function(data,indices){
  boot_data = data[indices, ]
  model4 <- glm(y~cig_cat, data=boot_data,family = "binomial")
  predicted4 <- predict(model4,boot_data, type="response")
  auc=as.numeric(pROC::auc(boot_data$y,predicted4,quiet=T))
  return(c(auc))
}
AUCSmokingPrsBoot = function(data,indices){
  boot_data = data[indices, ]
  model4 <- glm(y~cig_cat+prs, data=boot_data,family = "binomial")
  predicted4 <- predict(model4,boot_data, type="response")
  auc=as.numeric(pROC::auc(boot_data$y,predicted4,quiet=T))
  return(c(auc))
}

smokingmodel=function(allprs=GWAS24res$allprs)
{
  auc=data.frame(smoking=NA,smoking_low=NA,smoking_high=NA,prs=NA,prs_low=NA,prs_high=NA,
                 smokingprs=NA,smokingprs_low=NA,smokingprs_high=NA)
  
  OR1=data.frame(FORMER=NA,FORMER_low=NA,FORMER_high=NA,
                 OCCASIONAL=NA,OCCASIONAL_low=NA,OCCASIONAL_high=NA,
                 CURRENT=NA,CURRENT_low=NA,CURRENT_high=NA)
  OR2=data.frame(prs=NA,prs_low=NA,prs_high=NA)
  OR3=data.frame(prs=NA,prs_low=NA,prs_high=NA,
                 FORMER=NA,FORMER_low=NA,FORMER_high=NA,
                 OCCASIONAL=NA,OCCASIONAL_low=NA,OCCASIONAL_high=NA,
                 CURRENT=NA,CURRENT_low=NA,CURRENT_high=NA)
  pheno.prs=merge(phenotype[!is.na(phenotype$cig_cat),],allprs,by="ID")
  pheno.prs=pheno.prs[pheno.prs$ID %in% famtest$V2,]
  
  pheno.prs=myscale(pheno.prs)
  
  model1 <- glm(y~cig_cat, data=pheno.prs[pheno.prs$ID %in% famtest$V2 & !is.na(pheno.prs$cig_cat),],family = "binomial")
  coeff1 =coef(model1)
  OR1$FORMER=exp(coeff1)[2]
  OR1$OCCASIONAL=exp(coeff1)[3]
  OR1$CURRENT=exp(coeff1)[4]
  tmp1=confint(model1)
  OR1$FORMER_low=exp(tmp1[2,1])
  OR1$FORMER_high=exp(tmp1[2,2])
  OR1$OCCASIONAL_low=exp(tmp1[3,1])
  OR1$OCCASIONAL_high=exp(tmp1[3,2])
  OR1$CURRENT_low=exp(tmp1[4,1])
  OR1$CURRENT_high=exp(tmp1[4,2])
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$smoking=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCSmokingBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc$smoking_low=tmp$bca[4]
  auc$smoking_high=tmp$bca[5]
  
  model1 <- glm(y~prs, data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  coeff1 =coef(model1)
  OR2$prs=exp(coeff1)[2]
  tmp1=confint(model1)
  OR2$prs_low=exp(tmp1[2,1])
  OR2$prs_high=exp(tmp1[2,2])
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$prs=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc$prs_low=tmp$bca[4]
  auc$prs_high=tmp$bca[5]
  
  model1 <- glm(y~prs+cig_cat, data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  coeff1 =coef(model1)
  OR3$prs=exp(coeff1)[2]
  OR3$FORMER=exp(coeff1)[3]
  OR3$OCCASIONAL=exp(coeff1)[4]
  OR3$CURRENT=exp(coeff1)[5]
  tmp1=confint(model1)
  OR3$prs_low=exp(tmp1[2,1])
  OR3$prs_high=exp(tmp1[2,2])
  OR3$FORMER_low=exp(tmp1[3,1])
  OR3$FORMER_high=exp(tmp1[3,2])
  OR3$OCCASIONAL_low=exp(tmp1[4,1])
  OR3$OCCASIONAL_high=exp(tmp1[4,2])
  OR3$CURRENT_low=exp(tmp1[5,1])
  OR3$CURRENT_high=exp(tmp1[5,2])
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$smokingprs=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCSmokingPrsBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc$smokingprs_low=tmp$bca[4]
  auc$smokingprs_high=tmp$bca[5]
  
  # tmp=c(paste0("EV",1:10),"Age_cat","gender","prs","Age_cat")
  # idx=complete.cases(pheno.prs[,tmp])
  # table(idx)
  # tmp=paste0("y~",paste(c(paste0("EV",1:10),"Age_cat","gender","prs","Age_cat"),collapse="+"))
  # model1 <- glm(as.formula(tmp), data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  # predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  # auc$all=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  # 
  #model1 <- glm(y~prs+Age_cat+gender+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno.prs,family = "binomial")
  
  return(list(auc,OR1,OR2,OR3))
  
}
Smokingres=smokingmodel()

phenotype=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
phenotype$y=0
phenotype$y[which(phenotype$casecontrol=="CASE")]=1
idx=which(phenotype=="MISSING",arr.ind=T)
phenotype[idx]=NA
phenotype$cig_cat=factor(phenotype$cig_cat,levels=c("NEVER","FORMER","OCCASIONAL","CURRENT"))
methypheno=read_excel("../data/DNAmAge and Bladder Cancer.xlsx",sheet=4)
methypheno$Methylation=as.numeric(methypheno$`Methylation PC1`)
methypheno$cg05575921=as.numeric(methypheno$cg05575921)
methypheno$cg03636183=as.numeric(methypheno$cg03636183)
methypheno$cig_stat=gsub("Current Cigarette Smoker","Current",methypheno$cig_stat)
methypheno$cig_stat=gsub("Former Cigarette Smoker","Former",methypheno$cig_stat)
methypheno$cig_stat=gsub("Never Smoked Cigarettes","Never",methypheno$cig_stat)
methypheno$cig_stat=factor(methypheno$cig_stat,levels = c("Never","Former","Current"))
methypheno$tobacco_duration=as.numeric(methypheno$`tobacco duration (yrs)`)
methypheno$tobacco_cig_day=as.numeric(methypheno$`tobacco cig/day`)
methypheno$tobacco_cig_day[which(methypheno$cig_stat=="Never")]=0
methypheno$tobacco_duration[which(methypheno$cig_stat=="Never")]

tmp=phenotype[phenotype$study %in% c("ATBC","PLCO"),]
tmp=tmp[,c("study","study_pid","pid","id","gwas_id","sid","TGSID")]
write.csv(tmp,file="../result/PLCO_ATBC_sampleid.csv",row.names = F,quote=F)
methdat1=read.csv("../data/ImperialBladderCancer_PC1var.csv")
all(methdat1$TargetID==methypheno$TargetID)
cor(methdat1$PC1,methypheno$Methylation,use="complete")
methypheno$Methylation1=methdat1$PC1.28
table(is.na(methypheno$Methylation))
# FALSE  TRUE 
# 1598    40 
idx=which(!is.na(methypheno$Methylation))
table(methypheno$Study[idx],methypheno$CaseStatus[idx])
fm=glm(CaseStatus~cig_stat,data=methypheno,family = binomial)
coeff1 =coef(fm)
exp(coeff1)
fm=glm(CaseStatus~cig_stat,data=methypheno,family = binomial)
idx=methypheno$Study=="ATBC"
table(methypheno$cig_stat[idx],useNA = "ifany")
table(methypheno$cig_stat[!idx],useNA = "ifany")
table(methypheno$Study,methypheno$cig_stat)
table(methypheno$Study,methypheno$CaseStatus,useNA="ifany")
PLCOdat=methypheno[methypheno$Study=="PLCO",]
table(PLCOdat$cig_stat,PLCOdat$CaseStatus)
fm=glm(CaseStatus~cig_stat,data=methypheno[methypheno$Study=="PLCO",],family = binomial)
predicted1 <- predict(fm,methypheno, type="response")
auc1=as.numeric(pROC::auc(methypheno$CaseStatus,predicted1,quiet=T)) #0.524

coeff1 =coef(fm)
exp(coeff1)
fm=glm(CaseStatus~cig_stat+Study,data=methypheno,family = binomial)
coeff1 =coef(fm)
exp(coeff1)
tmp1=confint(fm)
exp(tmp1)

fm=glm("CaseStatus~Study",data=methypheno,family = binomial)
fm=glm("CaseStatus~Methylation",data=methypheno,family = binomial)
fm=glm("CaseStatus~Methylation1",data=methypheno,family = binomial)
summary(fm)$coefficients
# Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)  -0.0458885 0.05214911 -0.8799479 3.788875e-01
# Methylation1 -0.1395206 0.01665908 -8.3750497 5.520056e-17
predicted1 <- predict(fm,methypheno, type="response")
auc1=as.numeric(pROC::auc(methypheno$CaseStatus,predicted1,quiet=T))
fm=glm("CaseStatus~cg05575921",data=methypheno,family = binomial)
predicted2 <- predict(fm,methypheno, type="response")
auc2=as.numeric(pROC::auc(methypheno$CaseStatus,predicted2,quiet=T))
fm=glm("CaseStatus~cg03636183",data=methypheno,family = binomial)
predicted3 <- predict(fm,methypheno, type="response")
auc3=as.numeric(pROC::auc(methypheno$CaseStatus,predicted3,quiet=T))

fm=glm("CaseStatus~Methylation",data=PLCOdat,family = binomial)
fm=glm("CaseStatus~Methylation1",data=PLCOdat,family = binomial)
fm=glm("CaseStatus~cg05575921",data=PLCOdat,family = binomial)
fm=glm("CaseStatus~cg03636183",data=PLCOdat,family = binomial)
summary(fm)
cor(methypheno$cg05575921,methypheno$cg03636183,use="complete") #0.795

fm=glm(CaseStatus~cig_stat,data=methypheno,family = binomial)
fm=glm(CaseStatus~tobacco_duration,data=methypheno,family = binomial)
fm=glm(CaseStatus~tobacco_cig_day,data=methypheno,family = binomial)
fm=glm(CaseStatus~tobacco_duration+tobacco_cig_day,data=methypheno,family = binomial)
fm=glm(CaseStatus~cig_stat+tobacco_duration+tobacco_cig_day,data=methypheno,family = binomial)
predicted1 <- predict(fm,methypheno, type="response")
auc1=as.numeric(pROC::auc(methypheno$CaseStatus,predicted1,quiet=T))

fm=glm("CaseStatus~Methylation+Study",data=methypheno,family = binomial)
fm=glm("CaseStatus~cg05575921+Study+Age",data=methypheno,family = binomial)
fm=glm("CaseStatus~cg03636183+Study",data=methypheno,family = binomial)

cor(methypheno$cg05575921,methypheno$cg03636183,use="complete") #0.795
fm=glm("CaseStatus~cg05575921+cg03636183+Study+Age",data=methypheno,family = binomial)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.9096     0.3196   5.975  2.3e-09 ***
# cg05575921   -0.3835     0.5467  -0.701  0.48304    
# cg03636183   -2.9977     0.8932  -3.356  0.00079 ***
#Generate PRS for ATBC/PLCO samples
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
sumstatfile="../result/six10k_sumstats.txt"
CTselsnps=read.table("../result/CTselsnps.txt")$V1
sumdat=as.data.frame(fread(sumstatfile))
idx=match(CTselsnps,sumdat$rsid)
tmp=data.frame(snp=CTselsnps,effect=sumdat$a1[idx],beta=sumdat$beta[idx])
write.table(tmp,file="../result/CT.score",row.names = F,col.names = T,sep="\t",quote=F)
#compute PRS
cmd=paste0(plink2," --pfile ../result/imputation/merged_rsid  --score ../result/CT.score header cols=+scoresums --out ../result/CT_six10k")
system(cmd)
CTprs=as.data.frame(fread("../result/CT_six10k.sscore"))
table(CTprs$`#IID` %in% phenotype$study_pid)
idx=match(CTprs$`#IID`,phenotype$study_pid)
table(phenotype$study[idx])
mypheno=phenotype[phenotype$study %in% c("ATBC","PLCO"),]
mypheno=mypheno[mypheno$study_pid %in% CTprs$`#IID`,]
write.csv(mypheno[,c("study_pid","study")],file="../result/PLCO_ATBC_PRSID.csv",row.names = F)
write.csv(methypheno[,c("studyid","Study")],file="../result/PLCO_ATBC_methylationID.csv",row.names = F)

sum(mypheno$study_pid %in% methypheno$studyid)
smoking_methylation_model=function(allprs=GWAS24res$allprs)
{
  auc=data.frame(smoking=NA,smoking_low=NA,smoking_high=NA,prs=NA,prs_low=NA,prs_high=NA,
                 smokingprs=NA,smokingprs_low=NA,smokingprs_high=NA)
  
  OR1=data.frame(FORMER=NA,FORMER_low=NA,FORMER_high=NA,
                 OCCASIONAL=NA,OCCASIONAL_low=NA,OCCASIONAL_high=NA,
                 CURRENT=NA,CURRENT_low=NA,CURRENT_high=NA)
  OR2=data.frame(prs=NA,prs_low=NA,prs_high=NA)
  OR3=data.frame(prs=NA,prs_low=NA,prs_high=NA,
                 FORMER=NA,FORMER_low=NA,FORMER_high=NA,
                 OCCASIONAL=NA,OCCASIONAL_low=NA,OCCASIONAL_high=NA,
                 CURRENT=NA,CURRENT_low=NA,CURRENT_high=NA)
  pheno.prs=merge(phenotype,allprs,by="ID")
  pheno.prs=pheno.prs[pheno.prs$ID %in% famtest$V2,]
  pheno.prs=myscale(pheno.prs)
  
  model1 <- glm(y~cig_cat, data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  coeff1 =coef(model1)
  OR1$FORMER=exp(coeff1)[2]
  OR1$OCCASIONAL=exp(coeff1)[3]
  OR1$CURRENT=exp(coeff1)[4]
  tmp1=confint(model1)
  OR1$FORMER_low=exp(tmp1[2,1])
  OR1$FORMER_high=exp(tmp1[2,2])
  OR1$OCCASIONAL_low=exp(tmp1[3,1])
  OR1$OCCASIONAL_high=exp(tmp1[3,2])
  OR1$CURRENT_low=exp(tmp1[4,1])
  OR1$CURRENT_high=exp(tmp1[4,2])
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$smoking=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCSmokingBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc$smoking_low=tmp$bca[4]
  auc$smoking_high=tmp$bca[5]
  
  model1 <- glm(y~prs, data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  coeff1 =coef(model1)
  OR2$prs=exp(coeff1)[2]
  tmp1=confint(model1)
  OR2$prs_low=exp(tmp1[2,1])
  OR2$prs_high=exp(tmp1[2,2])
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$prs=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc$prs_low=tmp$bca[4]
  auc$prs_high=tmp$bca[5]
  
  model1 <- glm(y~prs+cig_cat, data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  coeff1 =coef(model1)
  OR3$prs=exp(coeff1)[2]
  OR3$FORMER=exp(coeff1)[3]
  OR3$OCCASIONAL=exp(coeff1)[4]
  OR3$CURRENT=exp(coeff1)[5]
  tmp1=confint(model1)
  OR3$prs_low=exp(tmp1[2,1])
  OR3$prs_high=exp(tmp1[2,2])
  OR3$FORMER_low=exp(tmp1[3,1])
  OR3$FORMER_high=exp(tmp1[3,2])
  OR3$CURRENT_low=exp(tmp1[4,1])
  OR3$CURRENT_high=exp(tmp1[4,2])
  predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  auc$smokingprs=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  boot_auc = boot(data =pheno.prs[pheno.prs$ID %in% famtest$V2,], statistic = AUCSmokingPrsBoot, R = 10000)
  tmp=boot.ci(boot_auc,type="bca")
  auc$smokingprs_low=tmp$bca[4]
  auc$smokingprs_high=tmp$bca[5]
  
  # tmp=c(paste0("EV",1:10),"Age_cat","gender","prs","Age_cat")
  # idx=complete.cases(pheno.prs[,tmp])
  # table(idx)
  # tmp=paste0("y~",paste(c(paste0("EV",1:10),"Age_cat","gender","prs","Age_cat"),collapse="+"))
  # model1 <- glm(as.formula(tmp), data=pheno.prs[pheno.prs$ID %in% famtest$V2,],family = "binomial")
  # predicted1 <- predict(model1,pheno.prs[pheno.prs$ID %in% famtest$V2,], type="response")
  # auc$all=as.numeric(pROC::auc(pheno.prs$y[pheno.prs$ID %in% famtest$V2],predicted1,quiet=T))
  # 
  #model1 <- glm(y~prs+Age_cat+gender+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, data=pheno.prs,family = "binomial")
  
  return(list(auc,OR1,OR2,OR3))
  
}

Smokingres=smoking_methylation_model()
save(CTres,CT5e8res,LDpredres,Lassosumres,PRScsres,GWAS24res,Smokingres,file="../result/PRS_res.RData")

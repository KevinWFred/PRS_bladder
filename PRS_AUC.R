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

phenotype=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
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
idx=match(pcadat$ID,phenotype$ID)
phenotype=cbind(phenotype,pcadat[idx,3:12])
idx=which(phenotype=="MISSING",arr.ind=T)
phenotype[idx]=NA
#CT

#p-value cutoffs (consistent with ../result/range_list)
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 

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
  sum(clumpedsnp$P<=pthres[idx_optimal]) #43 snps (P<5e-6) used
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
# [1] "auc:"
# tun       val
# 1 0.5691807 0.5954291
# auc_low  auc_high
# 1 0.5702395 0.5947593
# auc   auc_low  auc_high
# 1 0.5972794 0.5716296 0.5956955
#write.table(CTres$allprs,file="../result/CTprs.csv",sep=",",quote=F,row.names = F)

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


png("../result/CT_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
CTres_OR=prs_or_quantiles(allprs=CTres$allprs,main="CT",ytext = 0.27)
dev.off()
# or   or_low  or_high
# 1 1.40138 1.317046 1.492062

#LDpred2----------------------------
LDpred_prs=function(LDpredrda="../result/LDpred_six10k_pred.RData")
{
  load(LDpredrda)
  dim(df_beta)# 1004587 snps
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
# [1] "auc:"
# tun       val
# 1 0.5741286 0.5926333
# auc_low  auc_high
# 1 0.5714045 0.5949919
# auc   auc_low  auc_high
# 1 0.592982 0.5736006 0.5973283

#write.table(LDpredres$allprs,file="../result/LDpred_prs.csv",sep=",",quote=F,row.names = F)
png("../result/LDpred_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
LDpredres_OR=prs_or_quantiles(allprs=LDpredres$allprs,main="LDpred2",ytext=0.23)
dev.off()
# or   or_low or_high
# 1 1.395213 1.311029 1.48574

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
# [1] "auc:"
# tun       val
# 1 0.5694885 0.5913288
# auc_low  auc_high
# 1 0.56845 0.5923379
# auc   auc_low  auc_high
# 1 0.5915495 0.5693176 0.5934666

#write.table(Lassosumres$allprs,file="../result/Lassosum_prs.csv",sep=",",quote=F,row.names = F)
png("../result/Lassosum_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
Lassosumres_OR=prs_or_quantiles(allprs=Lassosumres$allprs,main="Lassosum2",ytext=0.2)
dev.off()
# or or_low or_high
# 1 1.395983 1.3114 1.48701

#PRS-CS
phi = c("phi1e+00","phi1e-02","phi1e-04","phi1e-06")
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

get_prscs=function(prefix="../result/prscs/six10k_pst_eff_a1_b0.5_",bimprefix="../result/six10k_train",outprefix="../result/six10k_train_prscs")
{
  if (!file.exists(paste0(outprefix,"PRS.sscore")))
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
  dim(allbeta) #1063497
  idx_optimal=which.max(auc_tun$EUR)
  allprs=data.frame(ID=c(PRS_tun$sample,PRS_val$sample),prs=c(PRS_tun[,idx_optimal+1],
               PRS_val[,idx_optimal+1]))
  aucs=get_aucs(allprs=allprs)
  return(list(aucs=aucs,allprs=allprs))
}
PRScsres=PRSCS_prs()
# auc
# tun       val
# 1 0.5583361 0.5827689
# 
# $auc_val
# auc_low  auc_high
# 1 0.5590809 0.5824027
# 
# $aucadj
# auc   auc_low  auc_high
# 1 0.5813707 0.5591336 0.5833984
png("../result/PRScs_prsquantiles.png",pointsize = 8,res=300,width=960,height=960)
PRScs_OR=prs_or_quantiles(allprs=PRScsres$allprs,main="PRS-CS",ytext=0.23)
dev.off()
# or   or_low  or_high
# 1 1.334052 1.254841 1.419003
write.table(PRScsres$allprs,file="../result/PRScs_prs.csv",sep=",",quote=F,row.names = F)

save(CTres,LDpredres,Lassosumres,PRScsres,file="../result/PRS_res.RData")

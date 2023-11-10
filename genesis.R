#!/usr/bin/env Rscript
#https://github.com/yandorazhang/GENESIS/blob/master/examples/examples.md

setwd("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/code")
.libPaths(c("/data/wangx53",.libPaths())) #where asian is installed

library(data.table)
library(GENESIS)
data(w_hm3.noMHC.snplist)

sumdat0=as.data.frame(fread("../result/six10k_sumstats.txt")) #n_eff=4/(1/case+1/contrl)
#used for genesis
sumdat=data.frame(snp=sumdat0$rsid,z=sumdat0$beta/sumdat0$beta_se,n=sumdat0$n_eff/4)
#QC
sumdat=preprocessing(sumdat)

fit2 <- genesis(sumdat[,1:3], filter=F, modelcomponents=2, cores=40, LDcutoff=0.1, LDwindow=1, c0=10, startingpic=0.005,qqplot.name="../result/twocom_qq.pdf")
fit2$estimates

est2 <- fit2$estimates$`Parameter (pic, sigmasq, a) estimates` # the model parameter estimates
v2 <- fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

starting <- rep(0,5)
starting[1] <- est2[1]
starting[2] <- 1/9
starting[3] <- est2[2]*5
starting[4] <- starting[3]/10
starting[5] <- est2[3]

fit3 <- genesis(sumdat[,1:3], filter=F, modelcomponents=3, cores=40, LDcutoff=0.1, LDwindow=1, c0=10,starting=starting,qqplot.name="../result/threecom_qq.pdf")
fit3$estimates

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

save(fit2,fit3,file="../result/GENESIS.RData")

#sbatch --mem=64g --cpus-per-task=41 --time=17:00:00 --gres=lscratch:64  /data/BB_Bioinformatics/Kevin/GWAS_Bladder/code/genesis.R

#estimate chip heritability
gwasdat2=as.data.frame(read_excel("../data/Supplemental Tables_R1_clean.xlsx",sheet=16,skip=2))
gwasdat2=gwasdat2[1:24,]
sumchip=0
for (i in 1:nrow(gwasdat2))
{
  se=log(gwasdat2$OR[i])/qnorm(gwasdat2$`P-value`[i]/2)
  sumchip=sumchip+2*gwasdat2$frequency[i]*(1-gwasdat2$frequency[i])*(log(gwasdat2$OR[i])^2-se^2)
}
sumchip #0.1908
fit3$estimates$`Total heritability in log-odds-ratio scale (sd)` #"0.63 (0.1283)"
sumchip/0.63 #0.303

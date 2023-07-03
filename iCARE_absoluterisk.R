#!/usr/bin/env Rscript
#absolute risk
library(iCARE)
library(dplyr)
library(ggplot2)
library(data.table)
famtrain=read.table("../result/six10k_train.fam")
famtest=read.table("../result/six10k_test.fam")

phenotype=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
changesampleid=function(id=tmp$`#IID`)
{
  dat=data.frame(id=id,study_pid=NA,TGSID=NA)
  tmp=intersect(dat$id,phenotype$study_pid)
  idx1=match(tmp,dat$id)
  idx2=match(tmp,phenotype$study_pid)
  dat$study_pid[idx1]=phenotype$study_pid[idx2]
  tmp=intersect(dat$id,phenotype$TGSID)
  idx1=match(tmp,dat$id)
  idx2=match(tmp,phenotype$TGSID)
  dat$study_pid[idx1]=phenotype$study_pid[idx2]
  if (sum(is.na(dat$study_pid))>0) stop("samples ID missing")
  return(dat$study_pid)
}
phenotype$ID=changesampleid(id=phenotype$study_pid)
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

#standardize based on controls
myscale=function(pheno.prs)
{
  idx=pheno.prs$casecontrol=="CONTROL"
  pheno.prs$prs=(pheno.prs$prs-mean(pheno.prs$prs[idx]))/sd(pheno.prs$prs[idx])
  return(pheno.prs)
}

opt="M"
load("../result/PRS_res.RData")
flcca_prs=CTres$allprs
pheno.prs.all=merge(phenotype,flcca_prs,by="ID")
pheno.prs.all=pheno.prs.all[pheno.prs.all$casecontrol=="CONTROL",]
if (opt=="M")
{
  pheno.prs=pheno.prs.all[pheno.prs.all$gender=="MALE",]
}else
{
  pheno.prs=pheno.prs.all[pheno.prs.all$gender=="FEMALE",]
}
pheno.prs=myscale(pheno.prs)
quantiles=quantile(pheno.prs$prs,c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1))
pheno.prs$prsquantile=cut(pheno.prs$prs,quantiles,labels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%"),include.lowest = T)
pheno.prs$prsquantile=factor(pheno.prs$prsquantile,levels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%"))

flcca_prs=data.frame(prs_sd=pheno.prs$prs,GWAS_ID=pheno.prs$ID,CASECONTROL_CASE=pheno.prs$casecontrol=="CASE",prs_quant=pheno.prs$prsquantile )

#incidence and mortality data prep

#load taiwan incidence rate file
#incidence_mortality_taiwan <- read.csv("../data/incidence_mortality_taiwan.csv")
# ages <- incidence_mortality_taiwan$Age
# incidence <- (incidence_mortality_taiwan$average_inc)/100000
# mortality <- (incidence_mortality_taiwan$mort_18)/100000

incidence_M =read.table("/data/BB_Bioinformatics/ProjectData/Bladder/Rates/Delay_Unadjusted/Male_Bladder_Cancer_Incidence_Rates.txt",sep="\t",header=T)
mortaility_M =read.table("/data/BB_Bioinformatics/ProjectData/Bladder/Rates/Delay_Unadjusted/Male_Bladder_Cancer_Mortality_Rates.txt",sep="\t",header=T)
incidence_F =read.table("/data/BB_Bioinformatics/ProjectData/Bladder/Rates/Delay_Unadjusted/Female_Bladder_Cancer_Incidence_Rates.txt",sep="\t",header=T)
mortaility_F =read.table("/data/BB_Bioinformatics/ProjectData/Bladder/Rates/Delay_Unadjusted/Female_Bladder_Cancer_Mortality_Rates.txt",sep="\t",header=T)
ages=50:84
if (opt=="M")
{
  incidence=as.numeric(incidence_M$Crude.Rate[incidence_M$Age.recode.with.single.ages.and.85.>=50 & incidence_M$Age.recode.with.single.ages.and.85.<=84])/100000
  mortality=as.numeric(mortaility_M$Crude.Rate[mortaility_M$Age.recode.with.single.ages.and.85.>=50 & mortaility_M$Age.recode.with.single.ages.and.85.<=84])/100000
}else
{
  incidence=as.numeric(incidence_F$Crude.Rate[incidence_F$Age.recode.with.single.ages.and.85.>=50 & incidence_F$Age.recode.with.single.ages.and.85.<=84])/100000
  mortality=as.numeric(mortaility_F$Crude.Rate[mortaility_F$Age.recode.with.single.ages.and.85.>=50 & mortaility_F$Age.recode.with.single.ages.and.85.<=84])/100000
}

lc_inc_gelac = cbind(ages, incidence)
lc_mort_gelac = cbind(ages, mortality)


#define iCARE parameters

prs_sd = list()
prs_sd[["name"]] <- "prs_sd"
prs_sd[["type"]] <- "continuous"

#define relative risk


lc_model_log_hr_prs = log(c(1.45)) #CT


lc_model_formula_prs = as.formula(diagnosis ~prs_sd)

lc_model_cov_info_prs = list(prs_sd)


#create vector of log (risk estimates)


lc_model_log_hr_prs <- setNames(lc_model_log_hr_prs, c("prs_sd"))


#create a data set with only controls to use as reference

gwas_con_ref <- flcca_prs %>%  
  filter(CASECONTROL_CASE==0) %>%
  select(c(prs_sd))
#N=2055(M) 714(F)

#run iCARE

set.seed(123)

alllevels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%")
computrisk=function(startage=50,interval=10,opt="start")
{
  #calculates lifetime risk from age 30 to 80
  
  res_dt = computeAbsoluteRisk(model.formula = lc_model_formula_prs,
                               model.cov.info = lc_model_cov_info_prs,
                               model.log.RR = lc_model_log_hr_prs,
                               model.ref.dataset = gwas_con_ref,
                               model.disease.incidence.rates = lc_inc_gelac, #matrix of age specific incidence rates
                               model.competing.incidence.rates = lc_mort_gelac, #matrix of age specific mortality rates
                               apply.age.start =startage ,
                               apply.age.interval.length = interval,
                               apply.cov.profile = gwas_con_ref,
                               return.refs.risk = TRUE)
  
  
  dt_risk <- as.data.frame(res_dt$details)
  dt_risk <- dt_risk %>%
    mutate(risk_p=Risk_Estimate *100)
  #all(dt_risk$prs_sd==gwas_con_ref) #T
  absrisk=cbind.data.frame(flcca_prs[flcca_prs$CASECONTROL_CASE==0,],dt_risk)
  #all(absrisk[,1]==absrisk[,8])#T
  
  endage=startage+interval
  if (opt=="start")
  {
    dat=data.frame(age=rep(startage,length(unique(flcca_prs$prs_quant))),mean_risk=NA,
                   prs_quant_flipped=alllevels)
  }else
  {
    dat=data.frame(age=rep(endage,length(unique(flcca_prs$prs_quant))),mean_risk=NA,
                   prs_quant_flipped=alllevels)
  }
  
  for (j in 1:length(alllevels))
  {
    idx=which(as.character(absrisk$prs_quant)==alllevels[j])
    dat$mean_risk[j]=mean(absrisk$Risk_Estimate[idx])
  }
  return(dat)
}

#10 years risk
risk10=NULL
for (i in 0:24)
{
  myrisk=computrisk(startage = 50+i)
  risk10=rbind(risk10,myrisk)
}
risk10$prs_quant_flipped=factor(risk10$prs_quant_flipped,levels = rev(alllevels))
final_dt_10year=risk10
#code for ggplot
final_dt_10year %>%
  ggplot(aes( x=age, y=mean_risk, group=as.factor(prs_quant_flipped), color=as.factor(prs_quant_flipped))) +
  geom_line() +
  xlab("Age (years)") + ylab("10-year Absolute Risk") +
  labs(color="")+
  
  #scale_color_discrete(labels=c(">95%", "90-95%", "80-90%", "70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "<5%")) +
  scale_color_discrete(labels=rev(alllevels)) +
  
  theme_bw() +
  theme(axis.text.x=element_text(colour="black",size=20),axis.text.y=element_text(colour="black",size=20),
        legend.text=element_text(size=10),
        legend.title = element_text(size=20),
        axis.title=element_text(size=20))
if (opt == "M")
{
  ggsave("../result/M_10yearrisk.pdf",width=8,height = 8)
}else
{
  ggsave("../result/F_10yearrisk.pdf",width=8,height = 8)
}


#cumulative risk
cumrisk=NULL
for (i in 51:84)
{
  myrisk=computrisk(startage = 50,interval = i-50,opt="end")
  cumrisk=rbind(cumrisk,myrisk)
}
cumrisk$prs_quant_flipped=factor(cumrisk$prs_quant_flipped,levels = rev(alllevels))

cumrisk %>%
  ggplot(aes( x=age, y=mean_risk, group=as.factor(prs_quant_flipped), color=as.factor(prs_quant_flipped))) +
  geom_line() +
  xlab("Age (years)") + ylab("Lifetime Absolute Risk") +
  labs(color="")+
  
  #scale_color_discrete(labels=c(">95%", "90-95%", "80-90%", "70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "<5%")) +
  scale_color_discrete(labels=rev(alllevels)) +
  
  theme_bw() +
  theme(axis.text.x=element_text(colour="black",size=20),axis.text.y=element_text(colour="black",size=20),
        legend.text=element_text(size=10),
        legend.title = element_text(size=20),
        axis.title=element_text(size=20))
if (opt == "M")
{
  ggsave("../result/M_cumrisk.pdf",width=8,height = 8)
}else
{
  ggsave("../result/F_cumrisk.pdf",width=8,height = 8)
}



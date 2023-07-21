#!/usr/bin/env Rscript
#absolute risk
.libPaths(c("/data/wangx53",.libPaths()))
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
idx=which(phenotype$cig_cat=="OCCASIONAL")
phenotype=phenotype[-idx,]
idx=which(phenotype$cig_cat=="")
phenotype=phenotype[-idx,]
idx=which(is.na(phenotype$cig_cat))
phenotype=phenotype[-idx,]
#phenotype$cig_cat=factor(phenotype$cig_cat,levels=c("NEVER","FORMER","CURRENT"))

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
#only consider validation samples
pheno.prs.all=pheno.prs.all[pheno.prs.all$ID %in% famtest$V2,]
#pheno.prs.all=pheno.prs.all[pheno.prs.all$casecontrol=="CONTROL",]
if (opt=="M")
{
  pheno.prs=pheno.prs.all[pheno.prs.all$gender=="MALE",]
}else
{
  pheno.prs=pheno.prs.all[pheno.prs.all$gender=="FEMALE",]
}

alllevels=c("<10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","> 90%")
pheno.prs=myscale(pheno.prs)
quantiles=quantile(pheno.prs$prs,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
pheno.prs$prsquantile=cut(pheno.prs$prs,quantiles,labels=alllevels,include.lowest = T)
pheno.prs$prsquantile=factor(pheno.prs$prsquantile,levels=alllevels)
if(sum(is.na(pheno.prs$prsquantile))>0) warning("some cases needs to process")
pheno.prs1 <- pheno.prs %>%
  group_by(prsquantile) %>%
  mutate(prs_mean=mean(prs)) %>%
  ungroup()

pheno.prs1$smoking_group=factor(pheno.prs1$cig_cat,levels=c("NEVER","FORMER","CURRENT"))
fm=glm("y~prs_mean+smoking_group+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10",data=pheno.prs1,family = binomial)
lc_model_log_hr_prs = summary(fm)$coefficients[c(3,4,2),1]#log(c(x,y,z))

#where x,y,z are the OR for former, current, and prs_mean, respectively.

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

#calculate prs_mean as the mean PRS value per each decile of PRS

#calculate mean PRS per decile

gelac_gwas_adeno <- flcca_prs %>%
  group_by(prs_quant) %>%
  mutate(prs_mean=mean(prs_sd)) %>%
  ungroup()
gelac_gwas_adeno <- as.data.frame(gelac_gwas_adeno)
flcca_prs=flcca_prs[!gelac_gwas_adeno$CASECONTROL_CASE,]
gelac_gwas_adeno=gelac_gwas_adeno[!gelac_gwas_adeno$CASECONTROL_CASE,]

idx=match(gelac_gwas_adeno$GWAS_ID,phenotype$study_pid)
gelac_gwas_adeno$smoking_group=phenotype$cig_cat[idx]

# Assuming variable “smoking_group” is categorized as the following
# 0=never; 1=former; 2=current
gelac_gwas_adeno$smoking_group[which(gelac_gwas_adeno$smoking_group=="NEVER")]=0
gelac_gwas_adeno$smoking_group[which(gelac_gwas_adeno$smoking_group=="FORMER")]=1
gelac_gwas_adeno$smoking_group[which(gelac_gwas_adeno$smoking_group=="CURRENT")]=2
gelac_gwas_adeno$smoking_group=factor(gelac_gwas_adeno$smoking_group,levels = c(0,1,2))

#define iCARE parameters
#code to add a factor variable to iCARE

smoking_group = list()
smoking_group[["names"]] <- "smoking_group"
smoking_group[["type"]] <- "factor"
smoking_group[["levels"]] <- c("0","1","2")


lc_model_formula_prs = as.formula(diagnosis ~ as.factor(smoking_group) +prs_mean)

prs_mean = list()
prs_mean[["name"]] <- "prs_mean"
prs_mean[["type"]] <- "continuous"
lc_model_cov_info_prs = list(smoking_group,prs_mean)


#create vector of log (risk estimates)


lc_model_log_hr_prs <- setNames(lc_model_log_hr_prs, c("smoking_group", "prs_mean"))

lc_model_log_hr_prs <- setNames(lc_model_log_hr_prs, c( "as.factor(smoking_group)1", "as.factor(smoking_group)2", "prs_mean"))

set.seed(123)


computrisk=function(startage=50,interval=30,opt="start")
{
  #calculates lifetime risk from age 30 to 80
  
  res_dt = computeAbsoluteRisk(model.formula = lc_model_formula_prs,
                               model.cov.info = lc_model_cov_info_prs,
                               model.log.RR = lc_model_log_hr_prs,
                               model.ref.dataset = gelac_gwas_adeno[,c(6,5)],
                               model.disease.incidence.rates = lc_inc_gelac, #matrix of age specific incidence rates
                               model.competing.incidence.rates = lc_mort_gelac, #matrix of age specific mortality rates
                               apply.age.start =startage ,
                               apply.age.interval.length = interval,
                               apply.cov.profile = gelac_gwas_adeno[,c(6,5)],
                               return.refs.risk = TRUE)
  
  
  dt_risk <- as.data.frame(res_dt$details)
  dt_risk <- dt_risk %>%
    mutate(risk_p=Risk_Estimate *100)
  #all(dt_risk$prs_sd==gwas_con_ref) #T
  absrisk=cbind.data.frame(flcca_prs,dt_risk)
  #all(absrisk[,1]==absrisk[,8])#T
  
 
  dat=data.frame(age=rep(startage,length(unique(flcca_prs$prs_quant))*3),smoking_group=rep(c(0,1,2),each=length(unique(flcca_prs$prs_quant))),mean_risk=NA,
                   prs_quant_flipped=alllevels)
  dat$smoking_group=factor(dat$smoking_group)
  
  
  for (j in 1:nrow(dat))
  {
    idx=which(as.character(absrisk$prs_quant)==dat$prs_quant_flipped[j] &absrisk$smoking_group==dat$smoking_group[j] )
    dat$mean_risk[j]=mean(absrisk$Risk_Estimate[idx])
  }
  colnames(dat)[ncol(dat)]="prs_cat10"
  return(dat)
}

gelac_gwas_risk_sum=computrisk()
gelac_gwas_risk_sum$cat10=as.numeric(factor(gelac_gwas_risk_sum$prs_cat10,levels = alllevels))
library(RColorBrewer)
#gelac_gwas_risk_sum$prs_cat10=as.numeric(gelac_gwas_risk_sum$prs_cat10)
ggplot(gelac_gwas_risk_sum, aes(x=cat10, y=mean_risk, fill=as.factor(smoking_group))) +
  geom_bar(stat = "identity",position=position_dodge())+
  ylab("Absolute risk (%)")+
  scale_x_continuous("Decile categories of PRS", labels = as.character(gelac_gwas_risk_sum$cat10), breaks = gelac_gwas_risk_sum$cat10)+
  theme_classic()+
  scale_fill_manual(name="", values = brewer.pal(n = 3, name = "Dark2"), label=c("Never smokers", "Former smokers","Current smokers")) +
  theme(axis.text.x=element_text(colour="black",size=20),axis.text.y=element_text(colour="black",size=20),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        axis.title=element_text(size=20),
        legend.position = "bottom")

if (opt == "M")
{
  ggsave("../result/M_cumrisk_smoking.pdf",width=10,height = 8)
}else
{
  ggsave("../result/F_cumrisk_smoking.pdf",width=10,height = 8)
}



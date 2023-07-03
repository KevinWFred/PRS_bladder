library(iCARE)
library(dplyr)
library(ggplot2)
#standardize based on controls
myscale=function(pheno.prs)
{
  idx=pheno.prs$CASECONTROL=="CONTROL"
  pheno.prs$prs=(pheno.prs$prs-mean(pheno.prs$prs[idx]))/sd(pheno.prs$prs[idx])
  return(pheno.prs)
}

famtrain=read.table("/data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EAS_never_train.fam")
famtest=read.table("/data/BB_Bioinformatics/Kevin/PRS_EASLC/result/EAS_never_test.fam")
load("/data/BB_Bioinformatics/Kevin/PRS_EASLC/result/PRS_fliccares1.RData")
allpres=CTSLEB_never$allprs
# phenotype2=phenotype[match(famtest$V1,phenotype$GWAS_ID),]
# pheno.prs=merge(phenotype2,flcca_prs,by="GWAS_ID")
flcca_prs=CTSLEB_never$allprs
pheno.prs=merge(phenotype,flcca_prs,by="GWAS_ID")
#pheno.prs$prs=scale(pheno.prs$prs)
pheno.prs=myscale(pheno.prs)
idx=which(pheno.prs$CASECONTROL=="CONTROL")
quantiles=quantile(pheno.prs$prs[idx],c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1))
pheno.prs$prsquantile=cut(pheno.prs$prs,quantiles,labels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%"),include.lowest = T)
pheno.prs$prsquantile=factor(pheno.prs$prsquantile,levels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%"))

#only for cases, not used
idx1=which(pheno.prs$prs>max(pheno.prs$prs[idx]))
if (length(idx1)>0)
{
  pheno.prs$prsquantile[idx1]=factor("> 95%",levels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%"))
}
idx1=which(pheno.prs$prs<min(pheno.prs$prs[idx]))
if (length(idx1)>0)
{
  pheno.prs$prsquantile[idx1]=factor("< 5%",levels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%"))
}

flcca_prs=data.frame(prs_sd=pheno.prs$prs,GWAS_ID=pheno.prs$GWAS_ID,CASECONTROL_CASE=pheno.prs$CASECONTROL=="CASE",prs_quant=pheno.prs$prsquantile )

#incidence and mortality data prep

#load taiwan incidence rate file
incidence_mortality_taiwan <- read.csv("../data/incidence_mortality_taiwan.csv")

ages <- incidence_mortality_taiwan$Age
incidence <- (incidence_mortality_taiwan$average_inc)/100000
mortality <- (incidence_mortality_taiwan$mort_18)/100000
lc_inc_gelac = cbind(ages, incidence)
lc_mort_gelac = cbind(ages, mortality)


#define iCARE parameters

prs_sd = list()
prs_sd[["name"]] <- "prs_sd"
prs_sd[["type"]] <- "continuous"

#define relative risk


lc_model_log_hr_prs = log(c(1.71)) #CT-SLEB


lc_model_formula_prs = as.formula(diagnosis ~prs_sd)

lc_model_cov_info_prs = list(prs_sd)


#create vector of log (risk estimates)


lc_model_log_hr_prs <- setNames(lc_model_log_hr_prs, c("prs_sd"))

lc_model_log_hr_prs <- setNames(lc_model_log_hr_prs, c( "prs_sd"))


#create a data set with only controls to use as reference

gwas_con_ref <- flcca_prs %>%  
  filter(CASECONTROL_CASE==0) %>%
  select(c( prs_sd))
#N=4544
#create a data set to be applied to 

#all
# gwas_all <- flcca_prs %>%  
#   select(c(prs_sd)) 


#create a data frame with GWAS id, case control status and quant category

# quant_dt <- flcca_prs %>%
#   dplyr::select(c(GWAS_ID, CASECONTROL_CASE, prs_quant)) %>%
#   filter(CASECONTROL_CASE==0)


#run iCARE

set.seed(123)

alllevels=c("< 5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-95%","> 95%")
computrisk=function(startage=30,interval=10,opt="start")
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
for (i in 0:40)
{
  myrisk=computrisk(startage = 30+i)
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
ggsave("../result/example_10yearrisk.pdf",width=8,height = 8)
#cumulative risk
cumrisk=NULL
for (i in 31:80)
{
  myrisk=computrisk(startage = 30,interval = i-30,opt="end")
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
ggsave("../result/example_cumrisk.pdf",width=8,height = 8)

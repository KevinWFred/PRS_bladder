#!/usr/bin/env Rscript

setwd("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/code")
.libPaths(c("/data/wangx53",.libPaths())) #where asian is installed

library(data.table)
library(GENESIS)

load("../result/GENESIS.RData")

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

projres=data.frame(n=c(seq(5000,200000,5000)),equal=NA,equal_low=NA,equal_high=NA,onem=NA,onem_low=NA,onem_high=NA,control10=NA,control10_low=NA,control10_high=NA)
sigma2toauc = function(x){ ifelse(x==0,0.50,pnorm(sqrt(0.5*x))) }
for (i in 1:nrow(projres))
{
  tmp=projection(est,v,n=projres$n[i]/2,CI=TRUE)
  projres$equal[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$equal_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$equal_high[i]=sigma2toauc(tmp$pheno.variance[3])
  #tmp=polyriskpredict(N=projres$n[i]/2, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=nrow(sumdat), M1=nrow(sumdat)*est[1], type="GWAS",alp.GWAS=5e-8, k.fold=3:5)
  # projres$equal[i]=tmp$AUC
  neff=1/(1/1e6 +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$onem[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$onem_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$onem_high[i]=sigma2toauc(tmp$pheno.variance[3])
  # tmp=polyriskpredict(N=neff, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=nrow(sumdat), M1=nrow(sumdat)*est[1], type="GWAS",alp.GWAS=5e-8, k.fold=3:5)
  # projres$onem[i]=tmp$AUC
  neff=1/(0.1/projres$n[i] +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$control10[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$control10_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$control10_high[i]=sigma2toauc(tmp$pheno.variance[3])
}
projct=projres

par(mar=c(4.5,5,2,1))
projres=data.frame(n=c(seq(5000,1000000,5000)),equal=NA,equal_low=NA,equal_high=NA,onem=NA,onem_low=NA,onem_high=NA,control10=NA,control10_low=NA,control10_high=NA,
                   control2=NA,control2_low=NA,control2_high=NA,control3=NA,control3_low=NA,control3_high=NA,control4=NA,control4_low=NA,control4_high=NA,
                   control5=NA,control5_low=NA,control5_high=NA)
sigma2toauc = function(x){ ifelse(x==0,0.50,pnorm(sqrt(0.5*x))) }
for (i in 1:nrow(projres))
{
  tmp=projection(est,v,n=projres$n[i]/2,CI=TRUE)
  projres$equal[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$equal_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$equal_high[i]=sigma2toauc(tmp$pheno.variance[3])
  #tmp=polyriskpredict(N=projres$n[i]/2, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=nrow(sumdat), M1=nrow(sumdat)*est[1], type="GWAS",alp.GWAS=5e-8, k.fold=3:5)
  # projres$equal[i]=tmp$AUC
  neff=1/(1/1e6 +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$onem[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$onem_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$onem_high[i]=sigma2toauc(tmp$pheno.variance[3])
  # tmp=polyriskpredict(N=neff, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=nrow(sumdat), M1=nrow(sumdat)*est[1], type="GWAS",alp.GWAS=5e-8, k.fold=3:5)
  # projres$onem[i]=tmp$AUC
  neff=1/(0.1/projres$n[i] +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$control10[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$control10_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$control10_high[i]=sigma2toauc(tmp$pheno.variance[3])
  neff=1/(0.5/projres$n[i] +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$control2[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$control2_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$control2_high[i]=sigma2toauc(tmp$pheno.variance[3])
  neff=1/(1/(projres$n[i]*3) +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$control3[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$control3_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$control3_high[i]=sigma2toauc(tmp$pheno.variance[3])
  neff=1/(1/(projres$n[i]*4) +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$control4[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$control4_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$control4_high[i]=sigma2toauc(tmp$pheno.variance[3])
  neff=1/(0.2/projres$n[i] +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$control5[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$control5_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$control5_high[i]=sigma2toauc(tmp$pheno.variance[3])
}
projCT=projres

projres=data.frame(n=c(seq(5000,400000,5000)),case30k=NA,case30k_low=NA,case30k_high=NA)
for (i in 1:nrow(projres))
{
  neff=1/(1/3e4 +1/projres$n[i])
  tmp=projection(est,v,n=neff,CI=TRUE)
  projres$case30k[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$case30k_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$case30k_high[i]=sigma2toauc(tmp$pheno.variance[3])
}
projCT30kcase=projres


#CT with observed data
projres=data.frame(n=c(res$n_eff*2,seq(100,10000,250)),equal=NA,equal_low=NA,equal_high=NA)
sigma2toauc = function(x){ ifelse(x==0,0.50,pnorm(sqrt(0.5*x))) }
for (i in 1:nrow(projres))
{
  tmp=projection(est,v,n=projres$n[i]/2,CI=TRUE)
  projres$equal[i]=sigma2toauc(tmp$pheno.variance[1])
  projres$equal_low[i]=sigma2toauc(tmp$pheno.variance[2])
  projres$equal_high[i]=sigma2toauc(tmp$pheno.variance[3])
  
}
projCTobs=projres

#res is AUC result based on meta analysis on combination of datasets, not ready
save(projCT,file="../result/projection.RData")
save(res,projCT,projCT30kcase,projct,projCTobs,projctobs,file="../result/projection.RData")
# 
# projres=data.frame(n=c(res$n_eff,seq(5000,100000,5000)),equal=NA,contrl10=NA,contrl1e6=NA)
# sigma2toauc = function(x){ ifelse(x==0,0.50,pnorm(sqrt(0.5*x))) }
# for (i in 1:nrow(projres))
# {
#   eq=projection(est,v,n=projres$n[i]/2,CI=TRUE)
#   contrl10=projection(est,v,n=1/(1/projres$n[i]+0.1/projres$n[i]),CI=TRUE)
#   contrl1e6=projection(est,v,n=1/(1/projres$n[i]+1e-6),CI=TRUE)
#   projres$equal[i]=sigma2toauc(eq$pheno.variance[1])
#   projres$contrl10[i]=sigma2toauc(contrl10$pheno.variance[1])
#   projres$contrl1e6[i]=sigma2toauc(contrl1e6$pheno.variance[1])
# }
# 
# plot(projres$n,projres$equal*1.04,type="l",lwd=3,ylab = "AUC",xlab="Case sample size in the GWAS for training PRS",col="red",cex.axis=1.2,cex.lab=1.2,pch=20,ylim=c(0.5,0.75))
# lines(projres$n,projres$contrl10*1.04,col="brown",lwd=3)
# lines(projres$n,projres$contrl1e6*1.04,col="blue",lwd=3)
# legend("bottomright",legend=c("Euqual Equal number of cases and controls", "Control case ratio 10:1","One Million controls"),col=c("red","brown","blue"),pch=16)


#draw plot

#Figures: Three curves LDPpred2 (1:1, 1:10, 1 million);  Three curves CT (1:1, 1:10, 1 million); Six figures separately with 95% CI; LDPred2 versus CT (given 1:1);
#1:1 case control ratio, CT and LDPred2 projection 95% CI with real data points

source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
library(ggplot2)
CIplot=function(data=projct,y="equal",high="equal_low",low="equal_high",title="CT, equal number of cases and controls",prefix="CT_CI_equl")
{
  scatplot=ggplot(data, aes(x = n, y = !!sym(y))) + geom_line(alpha = .5,color="#E64B35FF",linewidth=2)+
    labs(x = "Number of cases in the GWAS for training PRS",
         y = "AUC",
         linetype = "",
         title = title)+
    geom_ribbon(aes(ymin = !!sym(low), ymax = !!sym(high)), alpha = 0.1,fill="#E64B35FF") +
    theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_x_continuous(breaks=seq(5000,200000,5000))
  scatplot
  ggsave(filename=paste0("../result/",prefix,".png"),
         plot=scatplot, device="png",
         width=16, height=8, units="in", dpi=300)
}
CIplot(data=projct,title="CT, one million controls",prefix="CT_CI_onem")
CIplot(data=projct,y="onem",high="onem_low",low="onem_high",title="CT, one million controls",prefix="CT_CI_onem")
CIplot(data=projct,y="control10",high="control10_low",low="control10_high",title="CT, control case ratio 10:1",prefix="CT_CI_control10")

threeplot=function(data=projCT,title="CT",prefix="CT_3curve",ylim=c(0.5,0.75))
{
  scatplot=ggplot(data=data,aes(x=n))+geom_line(aes(y =equal,color="Equal number of cases and controls"), alpha=0.5,linewidth=2) +
    geom_line(aes(y = control4,color="Control case ratio 4:1"), alpha=0.5,linewidth=2) +
    geom_line(aes(y = onem,color="One million controls"), alpha=0.5,linewidth=2) +
    geom_line(aes(y = control10,color="Control case ratio 10:1"), alpha=0.5,linewidth=2) +
    scale_colour_manual("", 
                        breaks = c("Equal number of cases and controls", "Control case ratio 4:1","One million controls", "Control case ratio 10:1"),
                        values = c("#3C5488FF", "#00A087FF", "#F39B7FFF","#E64B35FF")) + ylim(ylim) +
    labs(x = "Number of cases in the GWAS for training PRS",
         y = "AUC",
         linetype = "",
         title = title)+
    theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_x_continuous(breaks=seq(5000,max(data$n),20000)) + theme(legend.position = "top")
  
  scatplot
  
  ggsave(filename=paste0("../result/",prefix,".png"),
         plot=scatplot, device="png",
         width=16, height=8, units="in", dpi=300)
}
threeplot(data=projCT,title="CT",prefix="CT_3curve")

threeplot1=function(data=projCT,title="CT",prefix="CT_3curve",ylim=c(0.5,0.75))
{
  scatplot=ggplot(data=data,aes(x=n))+geom_line(aes(y =equal,color="Equal number of cases and controls"), alpha=0.5,linewidth=2) +
    geom_line(aes(y = control3,color="Control case ratio 3:1"), alpha=0.5,linewidth=2) +
    geom_line(aes(y = onem,color="One million controls"), alpha=0.5,linewidth=2) +
    geom_line(aes(y = control10,color="Control case ratio 10:1"), alpha=0.5,linewidth=2) +
    scale_colour_manual("", 
                        breaks = c("Equal number of cases and controls", "Control case ratio 3:1","One million controls", "Control case ratio 10:1"),
                        values = c("#3C5488FF", "#00A087FF", "#F39B7FFF","#E64B35FF")) + ylim(ylim) +
    labs(x = "Number of cases in the GWAS for training PRS",
         y = "AUC",
         linetype = "",
         title = title)+
    theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_x_continuous(breaks=seq(5000,max(data$n),20000)) + theme(legend.position = "top")
  
  scatplot
  
  ggsave(filename=paste0("../result/",prefix,".png"),
         plot=scatplot, device="png",
         width=16, height=8, units="in", dpi=300)
}
threeplot1(data=projCT,title="CT",prefix="CT_3curve2")

equal2=data.frame(n=projct$n,ct=projct$equal,ct_low=projct$equal_low,ct_high=projct$equal_high,
                  ldpred=projld$equal,ldpred_low=projld$equal_low,ldpred_high=projld$equal_high)
plot2=function(data=equal2,prefix="CT_LDpred_CI")
{
  scatplot=ggplot(data,aes(x=n))+
    geom_line(aes(y =ct,color="CT"), alpha=0.5,linewidth=2)+
    geom_line(aes(y =ldpred,color="LDpred2"), alpha=0.5,linewidth=2)+
    #geom_ribbon(aes(ymin = ct_low, ymax = ct_high), alpha = 0.1,fill="#E64B35FF") +
    #geom_ribbon(aes(ymin = ldpred_low, ymax = ldpred_high), alpha = 0.1,fill="#00A087FF") +
    scale_colour_manual("", 
                        breaks = c("CT", "LDpred2"),
                        values = c("#E64B35FF", "#00A087FF")) +
    labs(x = "Number of cases in the GWAS for training PRS",
         y = "AUC",
         linetype = "",
         title = "Equal number of cases and controls")+
    theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_x_continuous(breaks=seq(5000,200000,5000)) + theme(legend.position = "top")
  scatplot
  ggsave(filename=paste0("../result/",prefix,".png"),
         plot=scatplot, device="png",
         width=16, height=8, units="in", dpi=300)
}
plot2()

equal2obs=data.frame(n=projctobs$n,ct=projctobs$equal,ct_low=projctobs$equal_low,ct_high=projctobs$equal_high,
                     ldpred=projldobs$equal,ldpred_low=projldobs$equal_low,ldpred_high=projldobs$equal_high)
#oonly keep 7 points
# idxrm=1:15
# idxrm=idxrm[!idxrm %in% idxkeep]
# equal2obs=equal2obs[-idxrm,]
# plot2obs=function(data=equal2obs,prefix="CT_LDpred_obs_CI")
# {
#   obsdat=res
#   obsdat$n=obsdat$n_eff*2
#   scatplot=ggplot(data,aes(x=n))+
#     geom_line(aes(y =ct,color="CT"), alpha=0.5,linewidth=2)+
#     geom_line(aes(y =ldpred,color="LDpred2"), alpha=0.5,linewidth=2)+
#     geom_ribbon(aes(ymin = ct_low, ymax = ct_high), alpha = 0.3,fill="#E64B35FF") +
#     geom_ribbon(aes(ymin = ldpred_low, ymax = ldpred_high), alpha = 0.3,fill="#00A087FF") +
#     geom_point(data=obsdat,aes(x=n,y=CTauc),size=2,color="#E64B35FF")+
#     geom_point(data=obsdat,aes(x=n,y=LDpredauc),size=2,color="#00A087FF")+
#     scale_colour_manual("",
#                         breaks = c("CT", "LDpred2"),
#                         values = c("#E64B35FF", "#00A087FF")) +
#     labs(x = "Number of cases in the GWAS for training PRS",
#          y = "AUC",
#          linetype = "",
#          title = "Equal number of cases and controls")+
#     theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     scale_x_continuous(breaks=seq(500,10000,500)) + theme(legend.position = "top")
#   scatplot
#   ggsave(filename=paste0("../result/",prefix,".png"),
#          plot=scatplot, device="png",
#          width=16, height=8, units="in", dpi=300)
# }
plotobs=function(data=equal2obs,prefix="CT_obs_CI",y="ct",ylow="ct_low",yhigh="ct_high",yobs="CTauc",title="CT",ylim=c(0.5,0.72))
{
  obsdat=res[idxkeep,] #7 points
  obsdat$n=obsdat$n_eff*2
  scatplot=ggplot(data,aes(x=n))+
    geom_line(aes(y =!!sym(y)),color="#E64B35FF", alpha=0.5,linewidth=2)+
    geom_ribbon(aes(ymin = !!sym(ylow), ymax = !!sym(yhigh)), alpha = 0.1,fill="#E64B35FF") +
    geom_point(data=obsdat,aes(x=n,y=!!sym(yobs)),size=2,color="#E64B35FF")+ ylim(ylim)+
    labs(x = "Number of cases in the GWAS for training PRS",
         y = "AUC",
         linetype = "",
         title = paste0(title,", equal number of cases and controls"))+
    theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_x_continuous(breaks=seq(500,10000,500)) + theme(legend.position = "top")
  scatplot
  ggsave(filename=paste0("../result/",prefix,".png"),
         plot=scatplot, device="png",
         width=16, height=8, units="in", dpi=300)
}
plotobs()
plotobs(data=equal2obs,prefix="LDpred2_obs_CI",y="ldpred",ylow="ldpred_low",yhigh="ldpred_high",yobs="LDpredauc",title="LDpred2")

#figure1, LDpred and CT, 1:1 and 10:1 ratio, color represents ratio and shape represents the method, no CIs, there will be four curves in total
fig1dat=data.frame(n=projct$n,ct=projct$equal,ct10=projct$control10,
                   ldpred=projld$equal,ldpred10=projld$control10)
scatplot=ggplot(fig1dat,aes(x=n))+
  geom_line(aes(y =ct,color="CT",linetype="1:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =ldpred,color="LDpred2",linetype="1:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =ct10,color="CT",linetype="10:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =ldpred10,color="LDpred2",linetype="10:1"), alpha=0.5,linewidth=2)+
  #geom_ribbon(aes(ymin = ct_low, ymax = ct_high), alpha = 0.1,fill="#E64B35FF") +
  #geom_ribbon(aes(ymin = ldpred_low, ymax = ldpred_high), alpha = 0.1,fill="#00A087FF") +
  scale_colour_manual("Method", 
                      breaks = c("LDpred2","CT"),
                      values = c("#E64B35FF", "#00A087FF")) +
  scale_linetype_manual("Control to case ratio",values=c("1:1"="solid","10:1"="longdash")) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 5, keyheight = 1),
         colour=guide_legend(keywidth = 3, keyheight = 1)) +
  labs(x = "Number of cases in the GWAS for training PRS",
       y = "AUC",
       linetype = "",
       title = "")+
  theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_continuous(breaks=seq(5000,200000,5000)) + theme(legend.position = "top")
scatplot
ggsave(filename=paste0("../result/projection_figure1_.png"),
       plot=scatplot, device="png",
       width=16, height=8, units="in", dpi=300)

#figure2 LDpred with 1:1, 2:1, 3:1, 4:1, 5:1, 10:1 ratio, color represents different ratio, add another red dashed vertical line at cases equal to 30K, no CIs, there are six curves and 1 dashed lines 3.
colfunc <- colorRampPalette(c("green", "blue"))
scatplot=ggplot(projld,aes(x=n))+
  geom_line(aes(y =equal,color="1:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =control2,color="2:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =control3,color="3:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =control4,color="4:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =control5,color="5:1"), alpha=0.5,linewidth=2)+
  geom_line(aes(y =control10,color="10:1"), alpha=0.5,linewidth=2)+
  geom_vline(xintercept = 30000, linetype="dotted", 
             color = "gray", size=1.5)+
  #geom_ribbon(aes(ymin = ct_low, ymax = ct_high), alpha = 0.1,fill="#E64B35FF") +
  #geom_ribbon(aes(ymin = ldpred_low, ymax = ldpred_high), alpha = 0.1,fill="#00A087FF") +
  scale_colour_manual("Control to case ratio", 
                      breaks = c("1:1", "2:1","3:1", "4:1", "5:1", "10:1"),
                      #values = c("#3C5488FF","#4DBBD5FF","#00A087FF","#8491B4FF","#F39B7FFF","#E64B35FF")) +
                      values=colfunc(6))+
  #scale_linetype_manual("Method",values=c("LDpred2"="solid","CT"="dashed")) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1,byrow=TRUE)) +
  labs(x = "Number of cases in the GWAS for training PRS",
       y = "AUC",
       linetype = "",
       title = "")+
  theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_continuous(breaks=seq(5000,200000,5000)) + theme(legend.position = "top")
scatplot
ggsave(filename=paste0("../result/projection_figure2_.png"),
       plot=scatplot, device="png",
       width=16, height=8, units="in", dpi=300)

#figure3 LDpred with cases fixed at 30K and different controls, with 95% CI
scatplot=ggplot(projld30kcase[projld30kcase$n<=200000,],aes(x=n))+
  geom_line(aes(y =case30k),color="#E64B35FF", alpha=0.5,linewidth=2)+
  geom_ribbon(aes(ymin = case30k_low, ymax = case30k_high), alpha = 0.1,fill="#E64B35FF") +
  geom_hline(yintercept = 0.66, linetype="dotted", 
             color = "gray", size=1.5)+
  scale_y_continuous(breaks=c(0.6,0.65,0.66,0.7))+
  labs(x = "Number of controls in the GWAS for training PRS",
       y = "AUC",
       linetype = "",
       title = "30,000 cases")+
  theme_Publication()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_continuous(breaks=seq(5000,200000,5000)) 
ggsave(filename=paste0("../result/projection_figure3.png"),
       plot=scatplot, device="png",
       width=16, height=8, units="in", dpi=300)
#!/usr/bin/env Rscript
library(data.table)
library("plotrix") #axis.break
library("RColorBrewer")
library("optparse")
library(readr)
library(dplyr)
library(ggplot2)
setwd("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/code")
#check sumdat
datafolder="/data/BB_Bioinformatics/ProjectData/Bladder/"
onem=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/1M/out.plink.rsq03.01.txt"))
length(unique(onem$CASES_Num)) #1
onem$CASES_Num[1] #1102
onem$CONTROLS_Num[1] #1043
two5m=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/2.5M/out.plink.rsq03.01.txt"))
length(unique(two5m$CASES_Num)) #1
two5m$CASES_Num[1] #723
two5m$CONTROLS_Num[1] #464
six10=as.data.frame(fread(paste0(datafolder,"Summaries/summary/Overall/610K/out.plink.rsq03.01.txt")))
length(unique(six10$CASES_Num)) #1
six10$CASES_Num[1] #3638
six10$CONTROLS_Num[1] #5419
six60w=as.data.frame(fread(paste0(datafolder,"Summaries/summary/Overall/660W/out.plink.rsq03.01.txt")))
length(unique(six60w$CASES_Num)) #1
six60w$CASES_Num[1] #1918
six60w$CONTROLS_Num[1] #1876
cnio=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/CNIO/Overall_a_ordered_matched.plink.rsq03.maf01.txt"))
colnames(cnio)
# [1] "CHR"          "BP"           "SNP"          "A1"           "A2"           "SE"           "OR"           "L95"         
# [9] "U95"          "P"            "GENO"         "Rsq"          "A1_EAF"       "CASES_Num"    "CONTROLS_Num"
length(unique(cnio$CASES_Num)) #1
cnio$CASES_Num[1] #1643
cnio$CONTROLS_Num[1] #695
omnix1=as.data.frame(fread(paste0(datafolder,"Summaries/summary/Overall/OmniX_1/out.plink.rsq03.01.txt")))
length(unique(omnix1$CASES_Num)) #1
omnix1$CASES_Num[1] #145
omnix1$CONTROLS_Num[1] #172
omnix2=as.data.frame(fread(paste0(datafolder,"Summaries/summary/Overall/OmniX_2/out.plink.rsq03.01.txt")))
length(unique(omnix2$CASES_Num)) #1
omnix2$CASES_Num[1] #1800
omnix2$CONTROLS_Num[1] #4939
onco=as.data.frame(fread(paste0(datafolder,"Summaries/summary/Overall/Oncoarray/out.plink.rsq03.01.txt")))
colnames(onco)
# [1] "CHR"          "BP"           "SNP"          "A1"           "A2"           "SE"           "OR"           "L95"         
# [9] "U95"          "P"            "INFO"         "Rsq"          "GENOTYPED"    "ALL_EAF"      "CASES_EAF"    "CONTROLS_EAF"
# [17] "ALL_Num"      "CASES_Num"    "CONTROLS_Num" "ID" 
length(unique(onco$CASES_Num)) #1
onco$CASES_Num[1] #399
onco$CONTROLS_Num[1] #403
decode=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/deCODE/DECODE.Bladder_Cancer_05052017.180516.matched.plink.rsq03.maf01.txt"))
colnames(decode)
# [1] "CHR"          "BP"           "SNP"          "A1"           "A2"           "SE"           "OR"           "L95"         
# [9] "U95"          "P"            "GENO"         "Rsq"          "A1_EAF"       "CASES_Num"    "CONTROLS_Num"
length(unique(decode$CASES_Num)) #1
decode$CASES_Num[1] #2079
decode$CONTROLS_Num[1] #327569

tmp=intersect(onco$SNP,decode$SNP)
idx1=match(tmp,onco$SNP)
idx2=match(tmp,decode$SNP)
cor(onco$ALL_EAF[idx1],decode$A1_EAF[idx2]) #0.99

#work on meta (without 660K)
#draw QQ plot
qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))       
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}

#optbreak=1, break top p-values, used for very low p-value cases
plotqq=function(data,optbreak=1,title="")
{
  dat = data %>% 
    mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>% 
    select(rsid,CHR,BP,P,MAF) %>% 
    rename(SNP = rsid)
  
  x = dat$P
  z = qnorm(x / 2)
  lambda = round(median(z^2) / qchisq(0.5,1), 3)
  N.effect = median(data$N)
  lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
  
  yLine <- c(-log10(5E-8))
  colLine <- c("red")
  dat$log10P = -log10(dat$P)
  gwas = as.data.frame(dat)
  # Determine frequency bins and create variable for binned QQ plot
  
  minMAF <- min(gwas$MAF)
  
  freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
  gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
  freqtable <- table(gwas$freqbin)
  freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
  freqtable <- freqtable[freqtable > 0]
  
  ## Generate QQ plot data by frequency bin
  fbin <- character(0)
  fN <- integer(0)
  fx <- numeric(0)
  fy <- numeric(0)
  fcol <- character(0)
  legendcol <- character(0)
  conf <- list()
  allcols <- brewer.pal(4,"Set1")
  ycol <- "log10P"
  for(f in 1:length(freqtable)){
    fbin <- c(fbin,names(freqtable)[f])
    fsnps <- which(gwas$freqbin ==names(freqtable)[f])
    plotdata <- qqplotdata(gwas[[ycol]][fsnps])
    fN <- c(fN,freqtable[f])
    fx <- c(fx,plotdata$e)
    fy <- c(fy,plotdata$o)
    fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
    conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                            'y'=c(plotdata$c975,rev(plotdata$c025)))
    legendcol <- c(legendcol,allcols[f])
  }
  legendtext <- paste0("MAF=",fbin,"; #SNPs=",format(fN,big.mark=",",scientific=FALSE))
  
  opt =  list(break.top = ifelse(optbreak==1,15,ceiling(max(fy))+1),
              top.size = 0.125)
  
  
  #png(filename = paste0(outpath,"/QQ_",eth[i1],"_",trait[i2],".png"), width = 8, height = 8, units = "in",res=300)
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)
  
  par(mar=c(5.1,5.1,4.1,1.1))
  
  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  #top
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  if (any(lab2>max(lab1))) lab2 <- lab2[lab2 > max(lab1)]
  
  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top
  
  # function to rescale the top part
  rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
  
  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
  # Plot confidence intervals	
  for(p in 1:length(conf)){
    polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
            col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
            border = NA)
  }
  
  # add points below top
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
  
  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  box()
  if (optbreak==1)
  {
    rescaled.y = rescale(fy[fy>opt$break.top])
    par(las=0)
    points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
    par(las=1)
    axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
    axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
    axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
    lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  }
  
  abline(h=ifelse(yLine<opt$break.top,
                  yLine,
                  rescale(yLine)),
         col=colLine,lwd=1.5,lty=2)
  legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
  text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
  text(5.9,1,paste(lambda_1000),cex = 1.5)
  
  title(title)
  #dev.off()
}
meta=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/result/Overall_no6101.tbl"))
#count sample size,add chr,pos, hg19
meta$ncontrols=meta$ncases=0
tmp=unlist(strsplit(meta$MarkerName,":"))
meta$CHR=as.integer(tmp[seq(1,length(tmp),2)])
meta$BP=as.integer(tmp[seq(2,length(tmp),2)])
meta$rsid=NA
updatemeta=function(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/1M/out.plink.rsq03.01.txt")
{
  dat=fread(myfile)
  snp=intersect(dat$SNP,meta$MarkerName)
  idx1=match(snp,dat$SNP)
  idx2=match(snp,meta$MarkerName)
  
  if(sum(colnames(dat) %in% c("ID"))>0)
  {
    meta$rsid[idx2]=dat$ID[idx1]
  }
  meta$ncontrols[idx2]=meta$ncontrols[idx2]+dat$CONTROLS_Num[idx1]
  meta$ncases[idx2]=meta$ncases[idx2]+dat$CASES_Num[idx1]
  
  return(meta)
}
meta=updatemeta()
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/2.5M/out.plink.rsq03.01.txt")
max(meta$ncontrols)
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/660W/out.plink.rsq03.01.txt")
max(meta$ncontrols)
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/CNIO/Overall_a_ordered_matched.plink.rsq03.maf01.txt")
max(meta$ncontrols)
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/OmniX_1/out.plink.rsq03.01.txt")
max(meta$ncontrols)
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/OmniX_2/out.plink.rsq03.01.txt")
max(meta$ncontrols)
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/Oncoarray/out.plink.rsq03.01.txt")
max(meta$ncontrols)
meta=updatemeta(myfile="/data/BB_Bioinformatics/ProjectData/Bladder/Summaries/summary/Overall/deCODE/DECODE.Bladder_Cancer_05052017.180516.matched.plink.rsq03.maf01.txt")
max(meta$ncontrols)
meta$N=1/(1/meta$ncases+1/meta$ncontrols)
write.table(meta,file="../result/metano610_sumstat.txt",sep="\t",quote=F,row.names = F)
meta1 = meta %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(Freq1),
         P = as.numeric(`P-value`),
         rsid = rsid,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
meta1=meta1[!is.na(meta1$rsid),]
meta1=meta1[meta1$FREQ_A1>0.01 & meta1$FREQ_A1<0.99,]
sum(meta1$P<5e-8) #394
png(filename = "../result/QQplot_metano610.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=meta1,optbreak=1,title="")
dev.off()

source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
plotmanhattan=function(data=meta1,title="",filename="../result/man_metano610.png")
{
  dat = data %>% 
    mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>% 
    select(rsid,CHR,BP,P,MAF) %>% 
    rename(SNP = rsid)
  dat$CHR1=factor(dat$CHR,levels = c(1:22,"X","Y"))
  idx=order(dat$CHR1,dat$BP)
  dat=dat[idx,]
  p.pwas <- 5E-08
  
  nCHR <- length(unique(dat$CHR))
  dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(dat$CHR)){
    nbp[i] <- max(dat[dat$CHR == i,]$BP)
    dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
    s <- s + nbp[i]
  }
  
  axis.set <- dat %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(dat$P)))) + 2 
  sig1 <- p.pwas
  
  
  sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
  library(ggplot2)
  manhplot <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                              color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.8, size=0.8) + 
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
    guides(color = FALSE) + 
    labs(x = NULL, 
         y = "-log10(p)", 
         linetype = "",
         title = title)+
    #subtitle = "A2: Critically ill COVID19+ vs. population controls;\nB1: Hospitalized COVID19+ vs non-hospitalized COVID19+;\nB2: Hospitalized COVID19+ vs. population controls;\nC2: Reported SARS-CoV-2 infection vs. population controls") + 
    theme_Publication()+
    theme(
      legend.position = "top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  
  #outpath <-"../result"  
  
  
  ggsave(filename=filename,
         plot=manhplot, device="png",
         width=9, height=4, units="in", dpi=300)
}
plotmanhattan()

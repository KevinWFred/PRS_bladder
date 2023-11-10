#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library(data.table)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
#work on 610K individual data
#tmp=fread("../result/imputation/chr22.traw",nrows=10000)
pheno=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
pheno$TGSID[which(pheno$TGSID=="#N/A")]=NA
#system("cp ../result/imputation/merged.psam ../result/imputation/merged.psam.orig")
#system("cp ../result/imputation/merged.fam ../result/imputation/merged.fam.orig")
#pick samples from 3 studies

six10dat=read.table("/data/BB_Bioinformatics/ProjectData/Bladder/Phenotype_data/Overall/610K/snptest.def",header=T)
six10dat=six10dat[-1,]
six10dat$ID=NA
idx=match(six10dat$ID_1,pheno$study_pid)
six10dat$ID=pheno$study_pid[idx]
idx1=which(is.na(idx))
idx2=match(six10dat$ID_1[idx1],pheno$TGSID)
six10dat$ID[idx1]=pheno$study_pid[idx2]
sum(is.na(six10dat$ID))

#restric to def file
idx=which(pheno$study %in% c("CPSII","NEBL","MDACC") & pheno$Illumina_Array=="610K")
table(pheno$casecontrol[idx])
all3studies=pheno$study_pid[idx]
all3studies=intersect(all3studies,six10dat$ID)
tmp=data.frame(FMID=0,IID=all3studies)
#write.table(tmp,file="../result/six10k_plinksample.txt",quote=F, row.names=F, col.names=F)

idx=match(all3studies,pheno$study_pid)
all3studies_case=all3studies[pheno$casecontrol[idx]=="CASE"]
all3studies_control=all3studies[!all3studies %in% all3studies_case]
set.seed(1000)
idx1=sample(1:length(all3studies_case),0.5*length(all3studies_case))
idx2=sample(1:length(all3studies_control),0.5*length(all3studies_control))
tunsamples=c(all3studies_case[idx1],all3studies_control[idx2])
# idx=sample(1:length(all3studies),0.5*length(all3studies))
# tunsamples=all3studies[idx]
valsamples=all3studies[!all3studies %in% tunsamples]
pheno$cig_cat=factor(pheno$cig_cat,levels=c("NEVER","FORMER","OCCASIONAL","CURRENT"))
pheno$y=0
pheno$y[which(pheno$casecontrol=="CASE")]=1
model1 <- glm(y~cig_cat, data=pheno[pheno$study_pid %in% valsamples,],family = "binomial")
summary(model1)
coeff1 =coef(model1)
exp(coeff1)
exp(confint(model1))


#transform to rsid
.libPaths(c("/data/wangx53",.libPaths()))
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
genome <- BSgenome.Hsapiens.UCSC.hg19
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37
seqlevelsStyle(genome) <- "NCBI"

tmp=as.data.frame(fread("../result/imputation/merged.bim"))
dat=data.frame(snp=tmp$V2,chr=tmp$V1,pos=tmp$V4)
findrsid=function(dat)
{
  dat$rsid=NA
  ## construct a GPos object containing all the positions we're interested in
  positions <- GPos(seqnames = dat$chr, pos = dat$pos)
  
  ## query the genome with out positions
  my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)
  
  ## this gives us a GPos object
  my_snps = as.data.frame(my_snps)
  idx = match(paste0(dat$chr,":",dat$pos),paste0(my_snps$seqnames,":",my_snps$pos))
  dat$rsid=my_snps$RefSNP_id[idx]
  dat$rsid=my_snps$RefSNP_id[idx]
  table(is.na(dat$rsid))
  return(dat)
}
dat=findrsid(dat)
tmp=data.frame(oldname=dat$snp,newname=dat$rsid)
tmp=tmp[!is.na(tmp$newname),]
write.table(tmp,file="../result/merged_updatename.txt",row.names = F,col.names = F,quote=F,sep="\t")
cmd=paste0(plink2," --pfile ../result/imputation/merged --update-name ../result/merged_updatename.txt --make-pgen --out ../result/imputation/merged_rsid")
system(cmd)

#generate training/testing data

trainsample=data.frame(FMID=0,IID=tunsamples)
trainidx=match(tunsamples,pheno$study_pid)
table(pheno$y[trainidx])
# 0    1 
# 1384 1307

testidx=match(valsamples,pheno$study_pid)
table(pheno$y[testidx])
# 0    1 
# 1385 1307
testsample=data.frame(FMID=0,IID=valsamples)
write.table(trainsample,file="../result/six10k_train_plinksample.txt",quote=F, row.names=F, col.names=F)
write.table(testsample,file="../result/six10k_test_plinksample.txt",quote=F, row.names=F, col.names=F)

# cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged_rsid --keep ../result/six10k_plinksample.txt --maf 0.01 ",
#           "--make-pgen --out ../result/six10k_maf01")
# system(cmd)

# cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/six10k_maf01 --freq ",
#           "--out ../result/six10k_maf01")
# system(cmd)
# tmp=fread("../result/six10k_maf01.afreq")
cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged_rsid --keep ../result/six10k_plinksample.txt --maf 0.01 ",
          "--make-bed --out ../result/six10k_maf01")
system(cmd)

cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged_rsid --keep ../result/six10k_train_plinksample.txt ",
          "--make-pgen --out ../result/six10k_train")
system(cmd)

cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged_rsid --keep ../result/six10k_test_plinksample.txt ",
          "--make-pgen --out ../result/six10k_test")
system(cmd)

cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --bfile ../result/six10k_maf01 --keep ../result/six10k_train_plinksample.txt --maf 0.01 --geno 0.05 ",
          "--make-bed --out ../result/six10k_train")
system(cmd)
tmp=fread("../result/six10k_train.bim")
write.table(tmp$V2,file="../result/six10k_train.snp",row.names = F,col.names = F,quote=F)
cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged_rsid --keep ../result/six10k_test_plinksample.txt ",
          "--extract ../result/six10k_train.snp --make-bed --out ../result/six10k_test")
system(cmd)
cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --bfile ../result/six10k_train --freq",
          " --out ../result/six10k_train")
system(cmd)
tmp=fread("../result/six10k_train.afreq")
quantile(tmp$ALT_FREQS)

tmp1=fread("../result/six10k_train.bim")
tmp2=fread("../result/six10k_test.bim")
all(tmp1$V5==tmp2$V5) #T
tmp3=fread("../result/six10k_test.pvar")
idx=match(tmp2$V2,tmp3$ID)
all(tmp2$V5==tmp3$ALT[idx]) #T
# tmp4=fread("../result/six10k_sumstats.txt")
# #to make the sumstats consistent with individual genotype (target data)
# idx=match(tmp2$V2,tmp4$rsid)
# all(tmp2$V5==tmp4$a1[idx],na.rm=T) #T

#to make sumstat for LDpred2
#six10k_sumstats.txt
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
formsumstats=function(sumstatfile="../result/metano610_sumstat.txt",
                      prefix_tun="../result/six10k_train",outprefix="six10k")
{
  sumdat=as.data.frame(fread(sumstatfile)) #7805414
  # dat=data.frame(snp=sumdat$MarkerName,chr=sumdat$CHR,pos=sumdat$BP)
  # dat=findrsid(dat)
  # idx=which(!grepl("^rs",sumdat$rsid))
  # sumdat$rsid[idx]=NA
  # table(sumdat$rsid==dat$rsid,useNA="ifany")
  # table(is.na(sumdat$rsid),is.na(dat$rsid))
  # sumdat$rsid=dat$rsid
  #idx=which(is.na(sumdat$rsid))
  #sumdat$rsid[idx]=paste0(sumdat$CHR[idx],sumdat$BP[idx])
  #sumdat=sumdat[!is.na(sumdat$rsid),]
  idx=which(colnames(sumdat) %in% c("chromosome","CHR"))
  if (length(idx)>0) colnames(sumdat)[idx]="chr"
  idx=which(colnames(sumdat) %in% c("position","BP"))
  if (length(idx)>0) colnames(sumdat)[idx]="pos"
  idx=which(colnames(sumdat) %in% c("SNP","rsid"))
  if (length(idx)>0) colnames(sumdat)[idx]="rsid"
  idx=which(colnames(sumdat) %in% c("alleleB","Effect.allele","Allele1"))
  if (length(idx)>0) colnames(sumdat)[idx]="a1" #effect
  idx=which(colnames(sumdat) %in% c("alleleA","Reference.allele","Allele2"))
  if (length(idx)>0) colnames(sumdat)[idx]="a0"
  
  sumdat$n_eff=4*sumdat$neff #For LDpred2
  
  idx=which(colnames(sumdat) %in% c("frequentist_add_beta_1","Effect"))
  if (length(idx)>0) colnames(sumdat)[idx]="beta"
  idx=which(colnames(sumdat) %in% c("frequentist_add_se_1","SE","StdErr"))
  if (length(idx)>0) colnames(sumdat)[idx]="beta_se"
  idx=which(colnames(sumdat) %in% c("frequentist_add_pvalue","P.value","P-value"))
  if (length(idx)>0) colnames(sumdat)[idx]="p"
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  sumstats=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=toupper(sumdat$a0),a1=toupper(sumdat$a1),n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  
  # # Open a temporary file
  # tmp <- tempfile(tmpdir = paste0("tmp-data",tmpdir))
  # on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  
  # preprocess the bed file (only need to do once for each data set)
  if (!file.exists(paste0(prefix_tun,".rds")))
    snp_readBed(paste0(prefix_tun,".bed"))
  # now attach the genotype object
  
  obj.bigSNP <- snp_attach(paste0(prefix_tun,".rds"))
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  #tmp=as.data.frame(fread(paste0(prefix_tun,".bim")))
  #tmp=fread("../data/FLCCA/plink2_EAS_never.pvar")
  # perform SNP matching
  info_snp <- snp_match(sumstats, map,join_by_pos=F)
  # idx=match(info_snp$rsid,sumstats$rsid)
  # sumstats1=sumstats[idx,]
  # table(sumstats1$beta==info_snp$beta)
  # table(sumstats1$a0==info_snp$a0)
  # all(sumstats1$beta_se==info_snp$beta_se)
  #7,805,414 variants to be matched.
  # 1,004,547 ambiguous SNPs have been removed.
  # 5,678,649 variants have been matched; 0 were flipped and 2,568,246 were reversed.
  sumstats=data.frame(chr=info_snp$chr,pos=info_snp$pos,rsid=info_snp$rsid,a0=info_snp$a0,a1=info_snp$a1,n_eff=info_snp$n_eff,beta_se=info_snp$beta_se,p=info_snp$p,beta=info_snp$beta)
  fwrite(sumstats,file=paste0("../result/",outprefix,"_sumstats.txt"),row.names = F,sep="\t",quote=F)
  #return(sumstats)
}

#24 GWAS snps
gwasdat1=as.data.frame(read_excel("../data/Supplemental Tables_R1_clean.xlsx",sheet=15,skip=3))
gwasdat2=as.data.frame(read_excel("../data/Supplemental Tables_R1_clean.xlsx",sheet=16,skip=2))
gwasdat2=gwasdat2[1:24,]
sum(gwasdat1$Position %in% gwasdat2$Position)
which(!gwasdat1$Position %in% gwasdat2$Position) #14
#gwasdat1[14,] #very close to another snps
idx=match(gwasdat2$Position,gwasdat1$Position)
gwasdat1=gwasdat1[idx,]
gwasdat2$SNP=gwasdat1$SNP
gwasdat=gwasdat2
#we explored a strong association signal for a low-quality imputed marker, rs36209093 (p = 3.21  1018; Supplementary Table 4, Supplementary Fig. 3). A proxy mar- ker (chr1:110229772) effectively tagged the GSTM1 deletion, improving the association signal (p = 8.84  1023)
#gwsdat$SNP[1] is not correct
#use ID instead
gwasdat$ID=paste0(gwasdat$chr,":",gwasdat$Position)
rm(gwasdat1,gwasdat2)

sumdat=as.data.frame(fread("../result/metano610_sumstat.txt"))
idx=which(!grepl("^rs",sumdat$rsid))
sumdat$rsid[idx]=paste0(sumdat$CHR[idx],":",sumdat$BP[idx])
idx=match(gwasdat$ID,sumdat$MarkerName)
quantile(sumdat$`P-value`[idx])
# 0%          25%          50%          75%         100% 
# 8.932000e-33 2.164125e-11 4.898000e-08 6.541250e-07 2.712000e-05

write.table(sumdat$rsid[idx],file="../result/Bladder_gwas24.snp",row.names = F,col.names=F,quote=F)
bim=as.data.frame(fread("../result/six10k_test.bim"))
which(bim$V2=="1:110229772") #1
tmp=data.frame(SNP=sumdat$rsid[idx],A1=sumdat$Allele1[idx],Effect=sumdat$Effect[idx])
tmp$A1=toupper(tmp$A1)
write.table(tmp,file="../result/metano610_score.txt",row.names = F,col.names = T,sep="\t",quote=F)
cmd=paste0("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged_rsid --extract ../result/Bladder_gwas24.snp --score ../result/metano610_score.txt cols=+scoresums,-scoreavgs header no-mean-imputation --out ../result/Bladder1_gwas24")
system(cmd)

#check on original 610K data
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
snppos=data.frame(chr=gwasdat$chr,start=gwasdat$Position,end=gwasdat$Position,name=gwasdat$SNP)
write.table(snppos,file="../result/Bladder_gwas24snp.range",row.names = F,col.names = F,sep=" ",quote=F)

extractsnp=function(infolder="/data/BB_Bioinformatics/ProjectData/Bladder/Imputed_data/610K/",
                    outfolder="../result/gwas_24snp/",prefix="Bladder_gwas24")
{
  allchrs=unique(snppos$chr)
  allchrs=allchrs[!allchrs %in% c("X","Y")]
  for (i in 1:length(allchrs))
  {
    chr=allchrs[i]
    cmd=paste0(plink2," --vcf ",infolder,"chr",chr,".dose.vcf.gz  --extract range ../result/Bladder_gwas24snp.range --memory 64000 --threads 8 --make-pgen --out ",outfolder,prefix,"_chr",chr)
    system(cmd)
  }
}
extractsnp()
mergedat=function(outfolder="../result/gwas_24snp/",prefix="Bladder_gwas24")
{
  allfiles=list.files(outfolder,paste0(prefix,"_chr\\w*.pvar"))
  tmp=data.frame(prefix=paste0(outfolder,allfiles))
  tmp$prefix=gsub(".pvar","",tmp$prefix)
  write.table(tmp,file=paste0(outfolder,prefix,"_merglist.txt"),row.names = F,col.names = F,quote=F)
  cmd=paste0(plink2," --pmerge-list ",outfolder,prefix,"_merglist.txt --make-pgen --out ",outfolder,prefix)
  system(cmd)
  #cmd=paste0(plink2," --pfile ",outfolder,prefix," --recode A-transpose --out ",outfolder,prefix)
  #system(cmd)
}#recover 22 snps
mergedat()
bim=read.table("../result/gwas_24snp/Bladder_gwas24.pvar")
idx=match(bim$V2,gwasdat$Position)
tmp=1:24
tmp[!tmp %in% idx] #0 are missing
tmp=data.frame(oldname=bim$V3,newname=gwasdat$SNP[idx])
tmp$newname[1]=tmp$oldname[1]
write.table(tmp,file="../result/gwas_24snp/Bladder_gwas24_updatename.txt",row.names = F,col.names = F,quote=F,sep="\t")
cmd=paste0(plink2," --pfile ../result/gwas_24snp/Bladder_gwas24 --update-name ../result/gwas_24snp/Bladder_gwas24_updatename.txt --make-pgen --out ../result/gwas_24snp/Bladder_gwas24")
system(cmd)
cmd=paste0("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/gwas_24snp/Bladder_gwas24 --extract ../result/Bladder_gwas24.snp --score ../result/metano610_score.txt cols=+scoresums,-scoreavgs header no-mean-imputation list-variants --out ../result/Bladder_gwas24")
system(cmd)
tmp=read.table("../result/Bladder_gwas24.sscore.vars")
tmp1=read.table("../result/Bladder_gwas24.sscore")
tmp2=read.table("../result/Bladder1_gwas24.sscore")
cor(tmp1$V4,tmp2$V5)
#[1] 0.9974624
prs_24gwas=read.table("../result/Bladder_gwas24.sscore")
cmd=paste0("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/gwas_24snp/Bladder_gwas24 --recode A-transpose --out ../result/gwas_24snp/Bladder_gwas24")
system(cmd)
runCT=function(sumstatfile="../result/six10k_sumstats.txt",prefix_tun="../result/six10k_train",prefix_val="../result/six10k_test",outprefix="six10k")
{
  #CT method
  #parameters for clumping
  pthr=1
  r2thr=0.1
  kbpthr=500
  
  cmd=paste0(plink," --bfile ",prefix_tun," --clump ",sumstatfile," --clump-p1 ",
             pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field rsid --clump-field p --out ../result/",outprefix)
  system(cmd)
  tmp=read.table(paste0("../result/",outprefix,".clumped"),header=T)
  write.table(tmp$SNP,file=paste0("../result/",outprefix,".clumpedsnp"),row.names=F,col.names = F,quote=F)
  sumstat=as.data.frame(fread(sumstatfile))
  tmp=data.frame(SNP=sumstat$rsid,A1=sumstat$a1,beta=sumstat$beta)
  write.table(tmp,file=paste0("../result/",outprefix,".score"),row.names=F,col.names=T,sep=" ",quote=F)
  tmp=data.frame(SNP=sumstat$rsid,P=sumstat$p)
  write.table(tmp,file=paste0("../result/",outprefix,".pvalue"),row.names=F,col.names=T,sep=" ",quote=F)
  cmd=paste0(plink2," --bfile ",prefix_tun," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ../result/",outprefix,"_tun")
  system(cmd)
  cmd=paste0(plink2," --bfile ",prefix_val," --score ../result/",outprefix,".score ",
             "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ../result/",outprefix,".pvalue --extract ../result/",outprefix,".clumpedsnp ",
             "--out ../result/",outprefix,"_val")
  system(cmd)
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 
  famtun=read.table(paste0(prefix_tun,".fam"))
  famval=read.table(paste0(prefix_val,".fam"))
  auc_tun=rep(0,length(pthres))
  for (i in 1:length(pthres))
  {
    prs=read.table(paste0("../result/",outprefix,"_tun.p_value_",i,".sscore"))
    if (any(!is.na(prs$V6)))
    {
      #all(prs$V1==famtun$V1)
      pheno.prs=data.frame(y=famtun$V6,prs=prs$V6)
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  
  idx_optimal=which.max(auc_tun)
  prs=read.table(paste0("../result/",outprefix,"_val.p_value_",idx_optimal,".sscore"))
  pheno.prs=data.frame(y=famval$V6,prs=prs$V6)
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  CTauc_val=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T)) #0.595
  return(CTauc_val)
}

tmp=runCT()

#LDpred2
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# 1. Read in the phenotype and covariate files
# read in phenotype and covariates
library(data.table)
library(magrittr)
run_ldpred2=function(tmpdir=1,sumstatfile="../result/six10k_sumstats.txt",prefix_tun="../result/six10k_train",prefix_val="../result/six10k_test",outprefix="six10k")
{
  print(sumstatfile)
  print(prefix_tun)
  print(prefix_val)
  #2. Obtain HapMap3 SNPs
  #LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs
  #load HapMap3 SNPs
  info <- readRDS(runonce::download_file(
    "https://ndownloader.figshare.com/files/25503788",
    fname = "map_hm3_ldpred2.rds"))
  
  #3. Load and transform the summary statistic file
  #Load summary statistic file
  # Read in the summary statistic file
  sumdat=as.data.frame(fread(sumstatfile)) #5678649
  print("sumdat dim:")
  print(dim(sumdat))
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  sumdat1=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=sumdat$a0,a1=sumdat$a1,n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  fwrite(sumdat1,file=paste0(sumstatfile,".ldpred"),row.names=F,sep="\t")
  sumstats <- bigreadr::fread2(paste0(sumstatfile,".ldpred")) 
  
  # Filter out hapmap SNPs
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,] #997853
  print("sumstat in HM3:")
  print(nrow(sumstats))
  #3. Calculate the LD matrix
  # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  tmp <- tempfile(tmpdir = paste0("tmp-data",tmpdir))
  #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  
  # preprocess the bed file (only need to do once for each data set)
  if (!file.exists(paste0(prefix_tun,".rds")))
    snp_readBed(paste0(prefix_tun,".bed"))
  # now attach the genotype object
  
  obj.bigSNP <- snp_attach(paste0(prefix_tun,".rds"))
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # perform SNP matching
  info_snp <- snp_match(sumstats, map)
  write.table(info_snp,file=paste0("../result/LDpred_",outprefix,"info_snp.txt"),row.names=F,sep="\t",quote=F)
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  genotype1 = snp_fastImputeSimple(genotype)
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  # calculate LD
  for (chr in 1:22) {
    print(chr)
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
      genotype,
      ind.col = ind.chr2,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      size = 3 / 1000
    )
    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp)
    } else {
      ld0=Matrix::colSums(corr0^2)
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
    if (sum(is.na(ld))>0) stop(chr)
  }
  save(ld,corr,file=paste0("../result/LDpred_",outprefix,"_ld.RData"))
  #4. Perform LD score regression
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  rownames(df_beta)=info_snp$rsid.ss
  ldsc <- snp_ldsc(   ld,
                      length(ld),
                      chi2 = (df_beta$beta / df_beta$beta_se)^2,
                      sample_size = df_beta$n_eff,
                      blocks = NULL)
  # ldsc <- snp_ldsc2(corr,df_beta)
  h2_est <- ldsc[["h2"]] #0.1
  print(paste0("h2_est:",h2_est))
  if (ldsc[['h2']] < 0) print('h2 negative')
  
  
  #6 grid model
  # Prepare data for grid model
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  h2_seq <- round(h2_est * c(0.3,0.7, 1, 1.4), 4)
  grid.param <-
    expand.grid(p = p_seq,
                h2 = h2_seq,
                sparse = c(FALSE, TRUE))
  # Get adjusted beta from grid model
  set.seed(1000) # to get the same result every time
  beta_grid <-
    snp_ldpred2_grid(corr, df_beta, grid.param, ncores = 1)
  
  
  famtun=read.table(paste0(prefix_tun,".fam"))
  pred_grid <- big_prodMat( genotype1, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)
  rownames(pred_grid)=famtun$V2
  
  #validation
  if(! file.exists(paste0(prefix_val,".rds")))
    snp_readBed(paste0(prefix_val,".bed"))
  val.obj.bigSNP <- snp_attach(paste0(prefix_val,".rds"))
  genotype_val <- val.obj.bigSNP$genotypes
  genotype1_val = snp_fastImputeSimple(genotype_val)
  famval=read.table(paste0(prefix_val,".fam"))
  pred_grid_val <- big_prodMat( genotype1_val, 
                                beta_grid, 
                                ind.col = info_snp$`_NUM_ID_`)
  rownames(pred_grid_val)=famval$V2
  
  save(ld,corr,beta_grid,pred_grid_val,beta_grid,pred_grid,ldsc,df_beta,grid.param,file=paste0("../result/LDpred_",outprefix,"_pred.RData"))
  
  #to get AUC
  famtun=read.table(paste0(prefix_tun,".fam"))
  famval=read.table(paste0(prefix_val,".fam"))
  auc_tun=rep(0,ncol(pred_grid))
  
  for (i in 1:ncol(pred_grid))
  {
    if (any(!is.na(pred_grid[,i])))
    {
      pheno.prs=cbind.data.frame(y=famtun$V6,prs=pred_grid[,i])
      
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  idx_optimal=which.max(auc_tun)
  param_optimal=grid.param[idx_optimal,]
  # 
  #auc on validation set
  
  pheno.prs=cbind.data.frame(y=famval$V6,prs=pred_grid_val[,idx_optimal])
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  LDpredauc_val=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))#0.593
  
  return(LDpredauc_val)
}

tmp1=run_ldpred2()
run_assosum2=function(tmpdir=1,sumstatfile="../result/six10k_sumstats.txt",prefix_tun="../result/six10k_train",prefix_val="../result/six10k_test",outprefix="six10k")
{
  print(sumstatfile)
  print(prefix_tun)
  print(prefix_val)
  #2. Obtain HapMap3 SNPs
  #LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs
  #load HapMap3 SNPs
  info <- readRDS(runonce::download_file(
    "https://ndownloader.figshare.com/files/25503788",
    fname = "map_hm3_ldpred2.rds"))
  
  #3. Load and transform the summary statistic file
  #Load summary statistic file
  # Read in the summary statistic file
  sumdat=as.data.frame(fread(sumstatfile)) #6081186
  print("sumdat dim:")
  print(dim(sumdat))
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  sumdat1=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=sumdat$a0,a1=sumdat$a1,n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  fwrite(sumdat1,file=paste0(sumstatfile,".ldpred"),row.names=F,sep="\t")
  sumstats <- bigreadr::fread2(paste0(sumstatfile,".ldpred")) 
  
  # Filter out hapmap SNPs
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,] #1005601
  print("sumstat in HM3:")
  print(nrow(sumstats))
  #3. Calculate the LD matrix
  # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  tmp <- tempfile(tmpdir = paste0("tmp-data",tmpdir))
  #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  
  # preprocess the bed file (only need to do once for each data set)
  if (!file.exists(paste0(prefix_tun,".rds")))
    snp_readBed(paste0(prefix_tun,".bed"))
  # now attach the genotype object
  
  obj.bigSNP <- snp_attach(paste0(prefix_tun,".rds"))
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # perform SNP matching
  info_snp <- snp_match(sumstats, map)
  write.table(info_snp,file=paste0("../result/LDpred_",outprefix,"info_snp.txt"),row.names=F,sep="\t",quote=F)
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  genotype1 = snp_fastImputeSimple(genotype)
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  # calculate LD
  for (chr in 1:22) {
    print(chr)
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
      genotype,
      ind.col = ind.chr2,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      size = 3 / 1000
    )
    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
    if (sum(is.na(ld))>0) stop(chr)
  }
  save(ld,corr,file=paste0("../result/LDpred_",outprefix,"_ld.RData"))
  #4. Perform LD score regression
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  #the above is the same as LDpred2
  rownames(df_beta)=info_snp$rsid.ss
  beta_lassosum2 <- snp_lassosum2(
    corr, df_beta, ncores = 1)
  
  famtun=read.table(paste0(prefix_tun,".fam"))
  pred_grid <- big_prodMat( genotype1, 
                            beta_lassosum2, 
                            ind.col = info_snp$`_NUM_ID_`)
  rownames(pred_grid)=famtun$V2
  
  #validation
  if(! file.exists(paste0(prefix_val,".rds")))
    snp_readBed(paste0(prefix_val,".bed"))
  val.obj.bigSNP <- snp_attach(paste0(prefix_val,".rds"))
  genotype_val <- val.obj.bigSNP$genotypes
  genotype1_val = snp_fastImputeSimple(genotype_val)
  famval=read.table(paste0(prefix_val,".fam"))
  
  pred_grid_val <- big_prodMat( genotype1_val, 
                                beta_lassosum2, 
                                ind.col = info_snp$`_NUM_ID_`)
  rownames(pred_grid_val)=famval$V2
  
  save(ld,corr,beta_lassosum2,pred_grid_val,pred_grid,df_beta,file=paste0("../result/Lassosum2_",outprefix,"_pred.RData"))
  
  #to get AUC
  famtun=read.table(paste0(prefix_tun,".fam"))
  famval=read.table(paste0(prefix_val,".fam"))
  auc_tun=rep(0,ncol(pred_grid))
  
  for (i in 1:ncol(pred_grid))
  {
    if (any(!is.na(pred_grid[,i])))
    {
      pheno.prs=cbind.data.frame(y=famtun$V6,prs=pred_grid[,i])
      
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  idx_optimal=which.max(auc_tun)
  
  # 
  #auc on validation set
  pheno.prs=cbind.data.frame(y=famval$V6,prs=pred_grid_val[,idx_optimal])
  model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  predicted1 <- predict(model1,pheno.prs, type="response")
  Lassosumauc_val=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))#0.591
  
  return(Lassosumauc_val)
}
tmp2=run_assosum2()


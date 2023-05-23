#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library(data.table)
#work on 610K individual data
#tmp=fread("../result/imputation/chr22.traw",nrows=10000)
pheno=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
#system("cp ../result/imputation/merged.psam ../result/imputation/merged.psam.orig")
#system("cp ../result/imputation/merged.fam ../result/imputation/merged.fam.orig")
tmp=fread("../result/imputation/merged.psam")
idx=match(tmp$`#IID`,pheno$TGSID)
idx1=which(is.na(idx))
idx2=match(tmp$`#IID`[idx1],pheno$study_pid)
idx[idx1]=idx2
sum(is.na(idx))
tmp$SEX=pheno$gender[idx]
tmp$SEX[which(tmp$SEX=="FEMALE")]=2
tmp$SEX[which(tmp$SEX=="MALE")]=1
tmp$Pheno=1
tmp$Pheno[which(pheno$casecontrol[idx]=="CASE")]=2
write.table(tmp,file="../result/imputation/merged.psam",sep=" ",quote=F,row.names = F)

tmp=fread("../result/imputation/merged.fam")
idx=match(tmp$V2,pheno$TGSID)
idx1=which(is.na(idx))
idx2=match(tmp$V2[idx1],pheno$study_pid)
idx[idx1]=idx2
sum(is.na(idx))
tmp$V5=pheno$gender[idx]
tmp$V5[which(tmp$V5=="FEMALE")]=2
tmp$V5[which(tmp$V5=="MALE")]=1
tmp$V6=1
tmp$V6[which(pheno$casecontrol[idx]=="CASE")]=2
write.table(tmp,file="../result/imputation/merged.fam",sep=" ",quote=F,row.names = F)

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
  table(is.na(dat$rsid))
  return(dat)
}
dat=findrsid(dat)
write.table(dat$snp[!is.na(dat$rsid)],file="../result/six60k_snpid.txt",row.names = F,col.names = F,quote=F)
cmd=paste0("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged --extract ../result/six60k_snpid.txt --make-pgen --out ../result/imputation/merged")
system(cmd)
tmp=as.data.frame(fread("../result/imputation/merged.pvar"))
idx=match(tmp$ID,dat$snp)
tmp$ID=dat$rsid[idx]
write.table(tmp,file="../result/imputation/merged.pvar", sep=" ",row.names = F,quote=F)

#generate training/testing data
tmp=fread("../result/imputation/merged.psam")
set.seed(1000) #sample0
trainidx=sample(1:nrow(tmp),floor(nrow(tmp)/2))
trainsample=data.frame(FMID=0,IID=tmp$`#IID`[trainidx])
table(tmp$Pheno[trainidx])
# 1    2 
# 2697 1831 

testidx=1:nrow(tmp)
testidx=testidx[!testidx %in% trainidx]
table(tmp$Pheno[testidx])
# 1    2 
# 2722 1807
testsample=data.frame(FMID=0,IID=tmp$`#IID`[testidx])
#EAS_never_train_plinksample.txt uses seed 1000
write.table(trainsample,file="../result/six10k_train_plinksample.txt",quote=F, row.names=F, col.names=F)
write.table(testsample,file="../result/six10k_test_plinksample.txt",quote=F, row.names=F, col.names=F)

cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --bfile ../result/imputation/merged --freq",
          "--out ../result/six10k")
system(cmd)
tmp=fread("../result/six10k.afreq")
cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged --keep ../result/six10k_train_plinksample.txt ",
          "--make-pgen --out ../result/six10k_train")
system(cmd)

cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../result/imputation/merged --keep ../result/six10k_test_plinksample.txt ",
          "--make-pgen --out ../result/six10k_test")
system(cmd)

cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --bfile ../result/imputation/merged --keep ../result/six10k_train_plinksample.txt ",
          "--make-bed --out ../result/six10k_train")
system(cmd)
cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --bfile ../result/imputation/merged --keep ../result/six10k_test_plinksample.txt ",
          "--make-bed --out ../result/six10k_test")
system(cmd)
cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --bfile ../result/six10k_train --freq",
          "--out ../result/six10k_train")
system(cmd)
tmp=fread("../result/six10k_train.afreq")
idx=which(tmp$ALT_FREQS<0.01 | tmp$ALT_FREQS>0.99)
tmp1=fread("../result/six10k_train.bim")
tmp2=fread("../result/six10k_test.bim")
all(tmp1$V5==tmp2$V5) #T
tmp3=fread("../result/six10k_test.pvar")
all(tmp1$V5==tmp3$ALT) #T

#to make the sumstats consistent with individual genotype (target data)

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
formsumstats=function(sumstatfile="../result/metano610_sumstat.txt",
                      prefix_tun="../result/six10k_train",outprefix="six10k")
{
  sumdat=as.data.frame(fread(sumstatfile)) #11404730
  dat=data.frame(snp=sumdat$MarkerName,chr=sumdat$CHR,pos=sumdat$BP)
  dat=findrsid(dat)
  table(sumdat$rsid==dat$rsid)
  sumdat$rsid=dat$rsid
  sumdat=sumdat[!is.na(sumdat$rsid),]
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
  idx=which(colnames(sumdat) %in% c("cases_total","N_cases","ncases"))
  if (length(idx)>0) colnames(sumdat)[idx]="ncases"
  idx=which(colnames(sumdat) %in% c("controls_total","N_controls","ncontrols"))
  if (length(idx)>0) colnames(sumdat)[idx]="ncontrols"
  sumdat$n_eff=4/(1/sumdat$ncases +1/sumdat$ncontrols)
  
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
  sumstats=data.frame(chr=info_snp$chr,pos=info_snp$pos,rsid=info_snp$rsid,a0=info_snp$a0,a1=info_snp$a1,n_eff=info_snp$n_eff,beta_se=info_snp$beta_se,p=info_snp$p,beta=info_snp$beta)
  write.table(sumstats,file=paste0("../result/",outprefix,"_sumstats.txt"),row.names = F,sep="\t",quote=F)
  #return(sumstats)
}

plink="/usr/local/apps/plink/1.9/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
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
    snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
  
  
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
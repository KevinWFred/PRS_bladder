#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library("sas7bdat")

dat=read.sas7bdat("../data/package-plco-957/Bladder File/plco_957_033022.sas7bdat")
table(dat$j_blad_cancer)
# 0      1 
# 152367   2520 
table(dat$j_blad_cancer_first)
# 0    1 
# 483 2037
table(dat$exc_blad_is_first_dx,dat$j_blad_cancer_first)
table(dat$exc_blad_is_first_dx)
# 0    1 
# 483 2037
table(dat$exc_in_tgwas_population)
# 0      1 
# 44325 110562 
table(dat$exc_in_tgwas_population,dat$exc_blad_is_first_dx)
#     0    1
# 0   99  569
# 1  384 1468

dat1=read.csv("../data/package-plco-957/Bladder File/plco_genotype_042022.csv")
table(dat1$population)
table(dat1$plco_id %in% dat$plco_id) 
# TRUE 
# 41261 
idx=match(dat1$plco_id,dat$plco_id)
table(dat1$age==dat$age[idx])
table(dat1$exc_blad_is_first_dx,dat$j_blad_cancer_first[idx])
pheno=read.csv("../data/PHENOTYPE_BLADDERGWAS_02102022_SK.csv")
pheno$TGSID[which(pheno$TGSID=="#N/A")]=NA
table(dat$plco_id %in% pheno$study_pid)

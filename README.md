# PRS for bladder cancer cancer in European population

- Goal  
The goal of the study is to build polygenic risk score (PRS) for bladder cancer in European population.  
- Data  
This is a large-scale study including GWAS summary data and/or individal genotype data from nine study/array groups, such as Spanish Bladder Cancer EPICURO Study (SBCS), European Prospective Investigation Into Cancer and Nutrition Study (EPIC),Alpha-Tocopherol, Beta-Carotene Cancer Prevention Study (ATBC),New England Bladder Cancer Study (NEBCS), Nijmegen Bladder Cancer Study (NBCS), Centro Nacional de Investigaciones Oncológicas (CNIO), DeCode Genetics, and others using individual genotypes from six genotyping platforms (llumina Human 1M Duo, Illumina Human 660W-Quad, Illumina Human 610-Quad, Illumina 2.5M, OncoArray, and Infinium OmniExpress).  
- Code  
read_gen.sh: read and QC on genotype data   
meta.R: meta analysis  
Run_PRS.R: the implement of C+T, LDpred2, and Lassosum2 methods  
run_prscs.py: the implement of PRS-CS method  
PRS_AUC.R: summarize AUC performance  

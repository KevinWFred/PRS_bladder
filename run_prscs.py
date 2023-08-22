#!/usr/bin/env python

import pandas as pd
import os
import numpy as np
import subprocess
import sys

os.chdir("/data/BB_Bioinformatics/Kevin/GWAS_Bladder/code")
os.getcwd()

#work on sumdat

#sumdatfile="../result/six10k_sumstats.txt"
#sumdat=pd.read_csv(sumdatfile,delimiter="\t")
#sumdat1=pd.DataFrame({"SNP":sumdat['rsid'],"A1":sumdat['a1'],"A2":sumdat['a0'],"BETA":sumdat["beta"],"P":sumdat['p']})
#sumdat1.to_csv("../result/six10k.sumdat.prscs",index=False,sep="\t")

print("start")
print(os.environ['OMP_NUM_THREADS'])
print(os.environ['HOSTNAME'])
print(os.environ['PWD'])
opt=int(sys.argv[1])
print(opt)
prscrsdir="/data/BB_Bioinformatics/Kevin/tools/PRScs/"
refdir="/data/BB_Bioinformatics/Kevin/tools/PRScsx/1kg/ldblk_1kg_eur"
#n_gwas from ../result/six10k_sumstats.txt, max(n_eff)
def run_prscs (ref_dir=refdir,bim_prefix="../result/six10k_train",
sst_file="../result/six10k.sumdat.prscs",
n_gwas="26632", OUTPUT_DIR="../result/prscs/six10k",
OUTPUT_FILE_PREFIX="six10k_train",
PARAM_PHI="1e-6",PHI="1e-06",SEED="1000"):
    #run on each chrom
    for chr in range(22):
        file1out=OUTPUT_DIR+OUTPUT_FILE_PREFIX+"_EUR_pst_eff_a1_b0.5_phi"+PHI+"_chr"+str(chr+1)+".txt"
        print(file1out)
        #if (not os.path.exists(file1out)):
        cmd=[prscrsdir+"PRScs.py", "--ref_dir="+ref_dir, "--bim_prefix="+bim_prefix,"--chrom="+str(chr+1), "--sst_file="+sst_file, "--n_gwas="+n_gwas, "--out_dir="+OUTPUT_DIR, "--phi="+PARAM_PHI, "--seed="+SEED]
        print(cmd)
        subprocess.call(cmd,shell=False)
    print("done")

if (opt==1):
    run_prscs()

if (opt==2):
    run_prscs(PARAM_PHI="1e-4",PHI="1e-04")

if (opt==3):
    run_prscs(PARAM_PHI="1e-2",PHI="1e-02")

if (opt==4):
    run_prscs(PARAM_PHI="1",PHI="1e+00")

tmp=pd.DataFrame({'code':['/data/BB_Bioinformatics/Kevin/GWAS_Bladder/code/run_prscs.py']*4,
                  'opt':[1,2,3,4]})
#tmp.to_csv("run_prscs.swarm",sep="\t",header=False,index=False)

#swarm --sbatch '--export=MKL_NUM_THREADS=1,NUMEXPR_NUM_THREADS=1,OMP_NUM_THREADS=1' -f /data/BB_Bioinformatics/Kevin/GWAS_Bladder/code/run_prscs.swarm --module python/3.9 -g 64 --time=0-15:00:00 --gres=lscratch:64 -p 2 

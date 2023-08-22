#!/usr/bin/env bash
#read 610k genotype data

plink=/usr/local/apps/plink/1.9/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2

cd /data/BB_Bioinformatics/Kevin/GWAS_Bladder/code
infolder="/data/BB_Bioinformatics/ProjectData/Bladder/Imputed_data/610K/"

$plink2 --vcf ${infolder}"chr22.dose.vcf.gz" dosage=DS --max-alleles 2 --make-pgen --out ../result/imputation/chr22 
$plink2 --pfile ../result/imputation/chr22 --make-bed --out ../result/imputation/chr22
$plink2 --bfile ../result/imputation/chr22 --recode A-transpose --out ../result/imputation/chr22

for chr in {1..22}
do
  echo $chr
  $plink2 --vcf ${infolder}chr${chr}".dose.vcf.gz" dosage=DS --max-alleles 2 --make-pgen --out ../result/imputation/chr$chr --threads 6 --memory 64000
  $plink2 --pfile ../result/imputation/chr$chr --make-bed --out ../result/imputation/chr$chr --threads 6 --memory 64000
done

outfolder="../result/imputation/"
mergechr(){
  local outfolder="$1"
  rm ${outfolder}mergelist.txt
  for chr in {1..22}
  do
    echo ${outfolder}chr${chr}  >> ${outfolder}mergelist.txt
  done
  
  $plink2 --pmerge-list ${outfolder}mergelist.txt --merge-max-allele-ct 2 --maf 0.005 --make-pgen --out ${outfolder}merged --memory 64000
}
$plink2 --pfile ../result/imputation/merged --maf 0.01 --make-bed --out ../result/imputation/merged
$plink2 --pfile ../result/imputation/merged --make-pgen --maf 0.01 --out ../result/imputation/merged

runPCA(){
  local bfile=$1
  local outfile=$2
  #$plink --bfile $bfile --maf 0.01 --geno 0.05  --make-bed --out ../result/$outfile

  # First, we need to perform prunning 
  $plink \
    --bfile $bfile \
    --indep-pairwise 200 50 0.25 \
    --out ../result/$outfile
  # Then we calculate the first 20 PCs
  $plink \
    --bfile $bfile \
    --extract ../result/$outfile.prune.in \
    --pca 20 \
    --threads 8\
    --out ../result/$outfile
}
bfile=../result/imputation/merged
outfile=merged

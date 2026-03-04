#!/bin/bash

#SBATCH -J gwas
#SBATCH -o e_o_files/male_add.%j.o
#SBATCH -e e_o_files/male_add.%j.e
#SBATCH -p skx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -A MCB25074
#SBATCH --mail-user=xliaoyi@my.utexas.edu
#SBATCH --mail-type=all

date=$1
pheno_path=$2
out_folder_name=$3
pheno_name=("${@:4}")  # Capture all arguments starting from the 4th one

plink2=../../EXECUTABLES/plink2


$plink2 \
    --glm no-x-sex hide-covar \
    --pfile ../../../ukb_pfile/wb_qced_maf0.001 \
    --extract ../../GWAS_INPUT/wb_qced_snps.snplist \
    --memory 256000 \
    --keep ../../GWAS_INPUT/20250703_both_WB_geno_qced_eids.txt \
    --keep-males \
    --pheno ${pheno_path} \
    --pheno-name ${pheno_name[@]} \
    --covar ../../GWAS_INPUT/20250719_covar_both.txt \
    --covar-name BIRTH_YEAR $(printf "PC%i " $(seq 1 10)) \
    --covar-variance-standardize \
    --out ../../PLINK_LM_OUTPUT/${out_folder_name}/${date}_male_add
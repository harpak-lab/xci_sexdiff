#!/bin/bash

#SBATCH -J magma_add
#SBATCH -o e_o_files/magma_add.%j.o
#SBATCH -e e_o_files/magma_add.%j.e
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=xliaoyi@my.utexas.edu
#SBATCH --mail-type=all

# binary traits sumstats contains extra "FIRTH?" column,
# so we need to have different code for q/b traits

date=$(date +%Y%m%d)

mkdir -p ../../MAGMA_OUTPUT/${date}_biomarkers

# qts=(
#   "HbA1c" "RBC_count" "albumin" "height" "bmi" "lymphocyte_perc"
# )

qts=(
  "testosterone" "SHBG" "glucose" "cholesterol" "vitamin_D" "IGF1"
  "total_protein" "oestradiol"
)

for qt in ${qts[@]}; do
  for sex in "male" "female"; do
    # sumstats=../../PLINK_LM_OUTPUT/all_13_traits_x/20250807_${sex}_X_add.${qt}.glm.linear
    sumstats=../../PLINK_LM_OUTPUT/20250923_biomarkers/20250923_${sex}_add.${qt}.glm.linear
    if [ -f ${sumstats} ]; then
      fn=$(basename ${sumstats} .glm.linear)
      echo ${fn}
      ../../EXECUTABLES/magma \
        --bfile ../../../ukb_bfile/ref_panel_10k \
        --gene-annot ../../MAGMA_INPUT/magma.genes.annot \
        --pval ${sumstats} use=3,15 ncol=11 \
        --genes-only \
        --out ../../MAGMA_OUTPUT/${date}_biomarkers/${fn}_window_5k
    fi
  done
done

bts=(
  "K50" "K51" "K80" "M05" "G43" "F20" "E050"
)

for bt in ${bts[@]}; do
  for sex in "male" "female"; do
    sumstats=../../PLINK_LM_OUTPUT/all_13_traits_x/20250807_${sex}_X_add.${bt}.glm.logistic.hybrid
    if [ -f ${sumstats} ]; then
      fn=$(basename ${sumstats} .glm.logistic.hybrid)
      echo ${fn}
      ../../EXECUTABLES/magma \
        --bfile ../../../ukb_bfile/ref_panel_10k \
        --gene-annot ../../MAGMA_INPUT/magma.genes.annot \
        --pval ${sumstats} use=3,16 ncol=12 \
        --genes-only \
        --out ../../MAGMA_OUTPUT/${date}_all_13_traits/${fn}_window_5k
    fi
  done
done

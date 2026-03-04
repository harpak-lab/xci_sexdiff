#!/bin/bash

# date=20251008
date=$(date +%Y%m%d)

pheno_name=(
  "HbA1c" "RBC_count" "albumin" "height" "bmi" "lymphocyte_perc"
  "K50" "K51" "K80" "M05" "G43" "F20" "E050"
)

pheno_path="../../GWAS_INPUT/20250911_phenos.txt"
out_folder_name="all_13_traits"

sbatch male_dom.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch male_add.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch female_dom.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch female_add.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch both_dom.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch both_add.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"

# for biomarker runs

date=$(date +%Y%m%d)
pheno_name=(
  "testosterone" "SHBG" "glucose" "cholesterol" "vitamin_D" "IGF1"
  "total_protein" "oestradiol"
)
pheno_path="../../GWAS_INPUT/20250923_phenos_biomarkers.txt"
out_folder_name="${date}_biomarkers"
mkdir -p ../../PLINK_LM_OUTPUT/${out_folder_name}

sbatch male_dom.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch male_add.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch female_dom.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch female_add.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch both_dom.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"
sbatch both_add.sh $date $pheno_path $out_folder_name "${pheno_name[@]}"

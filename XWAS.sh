#!/bin/sh

# phenotype command line argument
while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO 

# covariates include 10 PC, sex, and birth year; include in PHENO_DIR

# PARAMETERS AND PATHS
# paths and variables
QC_DIR=/QC
PHENO_DIR=/Phenotypes
OUTPUT_DIR=/Results

# beta = regression coefficient; orbeta = odds ratio
COL_NAMES=chrom,pos,ref,alt,omitted,test,nobs,beta,se,tz,p
COL_NAMES=chrom,pos,ref,alt,omitted,test,nobs,orbeta,se,tz,p


# PLINK2 GLM
# male
plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_DIR/ukb_imp_chrX_v3_11 --keep-males \
    --pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
    --covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 --covar-variance-standardize \
    --out $OUTPUT_DIR/male_X
    
# male dominant
plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar genotypic cols=$COL_NAMES \
    --pfile $QC_DIR/ukb_imp_chrX_v3_11 --keep-males \
    --pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
    --covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 --covar-variance-standardize \
    --out $OUTPUT_DIR/male_X_dom

#female
plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_DIR/ukb_imp_chrX_v3_11 --keep-females \
    --pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
    --covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 --covar-variance-standardize \
    --out $OUTPUT_DIR/female_X
    
# female dominant
plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar genotypic cols=$COL_NAMES \
    --pfile $QC_DIR/ukb_imp_chrX_v3_11 --keep-females \
    --pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
    --covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 --covar-variance-standardize \
    --out $OUTPUT_DIR/female_X_dom
    

    
    
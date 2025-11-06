
# make pgen
plink2 --memory 64000 --threads 16 --bgen $impute_path/ukb_imp_chrX_v3.bgen ref-first --sample $impute_path/ukb61666_imp_chrX_v3_s486605.sample \
    --make-pgen --out $qc_path/ukb_imp_chrX_v3_1
	
### GENOTYPE QC ###
# Info score >0.8
# obtain MFI ids from UKBB https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1967
plink2 --memory 64000 --threads 16 --pfile $qc_path/ukb_imp_chrX_v3_1 --extract $mfi_ids_path/ukb_mfi_chrX_v3_IDs.txt \
    --make-pgen --out $qc_path/ukb_imp_chrX_v3_2
	 
# call rate > 0.95 (missingness)
plink2 --memory 64000 --threads 16 --pfile $qc_path/ukb_imp_chrX_v3_2 --geno 0.05 --mind 0.05 --make-pgen --out $qc_path/ukb_imp_chrX_v3_3

# Alternate frequency (MAF) > 0.001 && <0.999
plink2 --memory 64000 --threads 16 --pfile $qc_path/ukb_imp_chrX_v3_3 --maf 0.001 --max-maf 0.999 --make-pgen --out $qc_path/ukb_imp_chrX_v3_4
    
# HWE > 1e-10
plink2 --memory 64000 --threads 16 --pfile $qc_path/ukb_imp_chrX_v3_4 --hwe 1e-10 --make-pgen --out $qc_path/ukb_imp_chrX_v3_5

# remove duplicates and keep only snps
plink2 --memory 64000 --threads 16 --pfile $qc_path/ukb_imp_chrX_v3_5 --snps-only --rm-dup exclude-all --make-pgen --out $qc_path/ukb_imp_chrX_v3_6

# remove indels and multiallelic snps
awk -F '\t' 'index($3, ":") !=0 {print $3}' $qc_path/ukb_imp_chrX_v3_6.pvar  > $qc_path/indels_chrX.txt
plink2 --pfile $qc_path/ukb_imp_chrX_v3_6 --exclude $qc_path/indels_chrX.txt --make-pgen --out $qc_path/ukb_imp_chrX_v3_7

### SAMPLE QC ###
# remove diff sex and sex chromosome aneuploidy samples
# UKBB genetic sex field ID: 22001
# UKBB reported sex field ID: 31
# UKBB aneuploidy field ID: 22019
plink2 --pfile $qc_path/ukb_imp_chrX_v3_7 --remove $pheno_path/diffsex_ids.txt $pheno_path/sex_aneuploidy_ids.txt --make-pgen --out $qc_path/ukb_imp_chrX_v3_8

# keep white british
# UKBB field ID: 22006
plink2 --pfile $qc_path/ukb_imp_chrX_v3_8 --keep $pheno_path/WBids.txt --make-pgen --out $qc_path/ukb_imp_chrX_v3_9

# keep unrelated individuals
# UKBB field ID: 22020
plink2 --pfile $qc_path/ukb_imp_chrX_v3_9 --keep $pheno_path/in_pca_ids.txt --make-pgen --out $qc_path/ukb_imp_chrX_v3_10

# remove withdrawn samples
plink2 --pfile $qc_path/ukb_imp_chrX_v3_10 --remove $pheno_path/withdrawn_ids.txt --make-pgen --out $qc_path/ukb_imp_chrX_v3_11



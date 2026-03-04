#!/bin/bash

#SBATCH -J magma_test
#SBATCH -o e_o_files/magma_test.%j.o
#SBATCH -e e_o_files/magma_test.%j.e
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=xliaoyi@my.utexas.edu
#SBATCH --mail-type=all

# annotate SNPs to genes

# you need to tell magma which SNP at which gene
# for --snp-loc, you can either provide a file which contains the SNP ID, chromosome, and position or a bim file in PLINK format
# for --gene-loc, you need to provide a gene location file you can download from https://cncr.nl/research/magma/ in Auxiliary files section

../../EXECUTABLES/magma \
    --annotate window=5,5 \
    --snp-loc ../../../ukb_bfile/wb_qced_maf0.001.bim \
    --gene-loc ../../MAGMA_INPUT/NCBI37.3.gene.loc \
    --out ../../MAGMA_INPUT/magma_window_5kb

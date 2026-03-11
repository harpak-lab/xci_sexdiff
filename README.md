# [Escape from X inactivation drives sex differences in gene expression](https://doi.org/10.1093/molbev/msag050)

Carrie Zhu<sup>1,2,+</sup>, Liaoyi Xu<sup>1,2,+</sup>, Arbel Harpak<sup>1,2,\*</sup>

<sup>1</sup> Department of Integrative Biology, University of Texas at Austin, Austin, TX, USA <br>
<sup>2</sup> Department of Population Health, University of Texas at Austin, Austin, TX, USA <br>
<sup>+</sup> C.Z. and L.X. contributed equally to this work <br>
<sup>*</sup> email: arbelharpak@utexas.edu 

----

Instructions to replicate data and figures in [Escape from X inactivation drives sex differences in gene expression](https://doi.org/10.1093/molbev/msag050) are provided below. Further explanation of process are detailed in Methods and Supplementary Text sections of paper.


## Software 
- R 4.3.1
- [plink 2.0](https://www.cog-genomics.org/plink/2.0/)


## ASE data
Take average ASE for genes with data from >1 individual from Gylemo et al. (2025) Table S4
- Input: elife-102701-supp1-v1.xlsx > convert to > Gylemo_ase_raw.txt
- Line 28-184 in fig1+4+s1-4+s6.Rmd
- Output: ASE_tissue_3_ind.txt

### Figure S5
Fisher exact test of individual x read count for 3 skewed XCI GTEx females, ASE data from Gylemo et al. (2025) Table S4
- Inputs: Gylemo_ase_raw.txt
- Run FigS5_Fisher.R


## GTEx Gene Expression
Use list of genes to get TPM data from GTEx v8
- Download GTEx v8 gct file: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct 
	- Bulk tissue expression: https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
- Gene lists from San Roman et. al.(1) Table S7
	- allNPX.txt
	- allPAR1.txt
- Run gene_filter_tacc.py 
- Output files
	- geneNPX_tpm_allsamples.txt
	- genePAR1_tpm_allsamples.txt

Filter and format data to include tissue and sex metadata, and remove male or female only tissues
- Download GTEx v8 sample attributes files
  - GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
	- GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
	- Metadata: https://www.gtexportal.org/home/downloads/adult-gtex/metadata
- Run format_TPM.R
	- Inputs: geneNPX_tpm_allsamples.txt and genePAR1_tpm_allsamples.txt
	- Outputs: geneNPX_TPM_filtered.txt and genePAR1_TPM_filtered.txt

Get female - male difference between median TPM per tissue per gene
- Run median_sex_TPM.R
	- Inputs: geneNPX_TPM_filtered.txt and genePAR1_TPM_filtered.txt
	- Outputs: meddiff_SMTSD_NPXgene.txt and meddiff_SMTSD_PAR1gene.txt

### Figure S1
SD ASE and Gene Expression Across Tissues
- Inputs: meddiff_SMTSD_NPXgene.txt and ASE_tissue_3_ind.txt
- Line 914-954 in fig1+4+s1-4+s6.Rmd 

### Figures 1; S2
Plot M-F median TPM difference against ASE
- Lines 183-344 (Fig 1) and 957-1001 (Fig S2) in fig1+4+s1-4+s6.Rmd
	- Inputs: meddiff_SMTSD_NPXgene.txt and meddiff_SMTSD_PAR1gene.txt and allNPX.txt


## XWAS
Quality Check 
- Run QC_X.sh
	- Inputs listed as datafields in script, requires UKBB access

Phenotypes
- Requires UKBB access
- Label as pheno_[pheno_name].txt
- 3 columns: #FID	IID	[pheno_name]

Perform XWAS using additive and dominance deviation models
- Run gwas.sh -p [pheno_name]
- XWAS files uploaded on https://www.harpaklab.com/data

### Figure 3
Examples of NPX SNPs with significant dominance effects
- Lines 27-44 in fig3+s7.ipynb

### Figure 4, S3, S4
- Gene Annotation
	- NCBI37.3.gene.loc in NCBI.zip
- Inputs
	- Xi expression data: XWAS summary statistics, ASE_tissue_3_ind.txt, NCBI37.3.gene.loc
Correlation between Xi expression and strength of association of trait (Fig S3)
- Line 1005-1128 in fig1+4+s1-4+s6.Rmd
Create heatmap of Xi expression and strength of association of trait across tissues (Fig 4)
- Line 347-709 in fig1+4+s1-4+s6.Rmd
Association between MAGMA p and xi expression heatmap (Fig S4)
- Line 712-913 in fig1+4+s1-4+s6.Rmd

### Figure S6 
Correlation between minimum gene-level p-value from trait association analysis and number of SNPs in gene
- Lines 1131-1457 in fig1+4+s1-4+s6.Rmd

### Figure S7
Comparison of the number of independent significant dominance-deviation SNPs across chromosomes
- "supp" in fig3+s7.ipynb

---

## References
1. Gylemo, B., et al. A whole-organism landscape of X-inactivation in humans. eLife 14, (2025). <https://doi.org/10.7554/eLife.102701.2>
3. San Roman, A. K. et al. The human inactive X chromosome modulates expression of the active X chromosome. Cell Genomics 3, 100259 (2023). <https://doi.org/10.1016/j.xgen.2023.100259>



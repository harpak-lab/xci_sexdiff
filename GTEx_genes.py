import pandas as pd

def get_genes(df_map, gene_list_input):
  # include only rows where GTEx file genes match those in gene_list_input
  # Description column is gene symbol
  return df_map[df_map['Description'].isin(gene_list_input)]

# load GTEx file
gtex_file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
df_full = pd.read_csv(gtex_file, sep="\t", usecols=["Description"])

# NPX
gene_list_npx = pd.read_csv("allNPX.txt", sep="\t", header=None)
gene_list_npx = gene_list_npx.iloc[:, 0].to_list()
df_npx = get_genes(df_full, gene_list_npx)
df_npx.to_csv("geneNPX_tpm_allsamples.txt", sep='\t', index=False)

# PAR1
gene_list_par1 = pd.read_csv("allPAR1.txt", sep="\t", header=None)
gene_list_par1 = gene_list_par1.iloc[:, 0].to_list()
df_par1 = get_genes(df_full, gene_list_par1)
df_par1.to_csv("genePAR1_tpm_allsamples.txt", sep='\t', index=False)

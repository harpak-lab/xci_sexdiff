import pandas as pd

df = pd.read_csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", usecols=["Description"])
gene_list = pd.read_csv("allNPX.txt", sep="\t", header=None)
gene_list = gene_list[0].to_list()
gene_id = df.index[df['Description'].isin(gene_list)].to_list()
gene_id = [x+1 for x in gene_id]
skipRows = list(range(1,56201))
skipRows = [x for x in skipRows if x not in gene_id]
df = pd.read_csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", skiprows=skipRows)
df.to_csv("geneNPX_tpm_allsamples.txt", sep='\t', index=False)

df = pd.read_csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", usecols=["Description"])
gene_list = ["AKAP17A", "ASMT", "ASMTL", "CD99", "CRLF2", "CSF2RA", "DHRSX", "GTPBP6", "IL3RA", "P2RY8", "PLCXD1", "PPP2R3B", "SHOX", "SLC25A6", "ZBED1"]
gene_list = gene_list[0].to_list()
gene_id = df.index[df['Description'].isin(gene_list)].to_list()
gene_id = [x+1 for x in gene_id]
skipRows = list(range(1,56201))
skipRows = [x for x in skipRows if x not in gene_id]
df = pd.read_csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", skiprows=skipRows)
df.to_csv("genePAR1_tpm_allsamples.txt", sep='\t', index=False)

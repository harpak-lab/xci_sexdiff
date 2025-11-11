#!/usr/bin/env Rscript

# REST API
library(httr)
library(jsonlite)
library(xml2)
library(tidyr)

fetch_endpoint <- function(geneId){
  print(geneId)
  server <- "https://rest.ensembl.org"
  ext <- paste0("/overlap/id/", geneId, "?feature=variation")
  r <- GET(paste(server, ext, sep = ""),  accept("application/json"))
  stop_for_status(r)
  sub_df <- as.data.frame(fromJSON(toJSON(content(r))))
  sub_df <- unnest(sub_df, cols=id)
  sub_df <- data.frame(id = sub_df$id, gene = rep(geneId, nrow(sub_df)))
  return(sub_df)
}

list_gene <- read.csv("gene_id_list.txt", sep="\t")

df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c("id", "gene")
for (i in 1:nrow(list_gene)) {
  geneId <- list_gene[i,1]
  sub_df <- fetch_endpoint(geneId)
  df <- rbind(df, sub_df)
}
head(list_gene)

gene_names <- read.csv("gene_list.txt", sep="\t")
df <- merge(gene_names, df, by.x = "ensembl_gene_id", by.y = "gene")

write.table(df, "ensemblGene_name_rsid.txt", sep="\t", row.names=FALSE, quote=FALSE)




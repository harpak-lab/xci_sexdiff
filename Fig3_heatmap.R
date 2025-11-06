library(ggplot2)
library(ggrepel)
library(reshape2)
library(stringr)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyr)
library(ggdendro)
library(Hmisc)
library(ggnewscale)

## load and format correlation data
# -logP and ASE correlation by tissue and phenotype
agg_df <- read.csv("allPheno_ase_bytissue.txt", sep="\t")

# p-value labels
agg_df <- agg_df %>%
  mutate(p_fdr = p.adjust(cor_pvalue, "fdr")) %>%
  mutate(cor_pvalue = case_when(cor_pvalue > 0.05 ~ NaN,    
                                TRUE ~ cor_pvalue)) %>%
  mutate(p_label = ifelse(!is.na(cor_pvalue), sapply(cor_pvalue, format, digits=1, nsmall=3), ""))

# abbreviate long names
agg_df <- agg_df %>%
  mutate(Name = case_when(Name=="Schizophrenia, schizotypal and delusional disorders" ~ "Schizophrenia",
                          Name=="Inflammatory polyarthropathies" ~ "Inflammatory \n polyarthropathies",
                          Name=="Thyrotoxicosis w/ diffuse goiter" ~ "Thyrotoxicosis",
                          TRUE ~ Name)) %>%
  mutate(Tissue = case_when(Tissue=="Skin - Not Sun Exposed (Suprapubic)" ~ "Skin - Suprapubic",
                            Tissue=="Cells - EBV - transformed lymphocytes" ~ "LCL",
                            TRUE ~ Tissue))

## Dendrogram
# pairwise correlation of ASE for tissues 
df_ase <- read.csv("ASE_tissue.txt", sep="\t")

df_ase <- df_ase %>%
  select(Gene, Tissue, ASE) %>%
  arrange(Tissue) %>%
  pivot_wider(names_from = Tissue, values_from = ASE) %>%
  as.data.frame()
row.names(df_ase) <- df_ase[,1]
df_ase[,1] <- NULL

# correlation matrix 
corr_matrix <- rcorr(as.matrix(df_ase), type="pearson")
tissue_dendro <- as.dendrogram(hclust(d = dist(x = corr_matrix$r)))

# dendrogram plot
dendro_plot <- ggdendrogram(data = tissue_dendro, rotate=T, labels = F, leaf_labels = F)
dendro_order <- order.dendrogram(tissue_dendro)

## Split to add and dom 
# split dataframe
separate_df <- function(mode) {
  df_half <- agg_df %>%
    filter(TEST == mode) %>%
    select(Tissue, cor_pvalue, p_fdr, Name, p_label, correlation) %>%
    mutate(correlation = ifelse(is.na(cor_pvalue), NaN, correlation)) %>%
    mutate(correlation = ifelse(is.na(p_fdr), NaN, correlation)) %>%
    as_tibble()
  df_half$Tissue <- factor(df_half$Tissue, levels=df_half$Tissue[dendro_order], ordered=T)
  return(df_half)
}
# create triangle polygon for plot
triangle_df <- function(df_half, mode) {
  m <- ifelse(mode=="ADD", 1, 0)
  df_half_poly <- df_half[rep(seq(nrow(dom_df)), each = 3),]
  df_half_poly <- df_half_poly %>%
    mutate(Name1 = as.numeric(as.factor(df_half_poly$Name))) %>%
    mutate(Tissue1 = as.numeric(as.factor(df_half_poly$Tissue))) %>%
    mutate(Name1 = Name1 + c(0.5-m, 0.5, -0.5)) %>%
    mutate(Tissue1 = Tissue1 + c(-0.5, 0.5-m, 0.5)) %>%
    mutate(z = rep(seq(nrow(df_half_poly)/3), each = 3))
  return(df_half_poly)
}
add_df <- separate_df("ADD"); dom_df <- separate_df("DOMDEV")
add_poly <- triangle_df(add_df, "ADD"); dom_poly <- triangle_df(dom_df, "DOMDEV")

add_df <- as.data.frame(add_df)
write.table(add_df, "fig3_data.txt", sep="\t", row.names = F, quote = F)

## PLOT HEATMAP
# heatmap
hm_plot <- ggplot(add_df, aes(x=Name, y=as.factor(Tissue))) +
  geom_tile(fill="gray90") +
  geom_polygon(data = add_poly, aes(x=Name1, y=Tissue1, group=z, fill=cor_pvalue), linewidth=0.2, color="gray90") +      # additive lower left triangle
  geom_text(aes(label=p_label), size=2.2, hjust=0.8, vjust=1.5) +
  geom_polygon(data = dom_poly, aes(x=Name1, y=Tissue1, group=z, fill=cor_pvalue), linewidth=0.2, color="gray90") +      # dominance upper right triangle
  geom_text(data = dom_df, aes(label=p_label), size=2.2, hjust=0.2, vjust=-0.5) +
  scale_fill_gradientn(colors = c("steelblue", "#ecf2f7"), na.value="white", name="Pearson p-value",
                       limits = c(0,0.05)) +
  scale_x_discrete(position="bottom") +
  labs(title = "Pearson correlation between ASE and XWAS -log10(p-values)") +
  theme(legend.position = "top", axis.title = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1, size=8, vjust=1.03), axis.text.y = element_text(size=8)) +
  coord_fixed(ratio=4/7)
hm_plot
grid.newpage()
print(hm_plot,
      vp = viewport(x = 0.42, y = 0.5, width = 0.85, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.90, y = 0.43, width = 0.15, height = 0.72))




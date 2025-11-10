library(ggplot2)
library(ggrepel)
library(reshape2)
library(stringr)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyr)

### FUNCTIONS AND DATA ###

# ASE and associated SNPs
df_ase <- read.csv("ASE_tissue.txt", sep="\t") 
colnames(df_ase) <- c("Gene", "Tissue", "ASE")

# get nearby associated genes
gene_list <- read.csv("ensemblGene_name_rsid.txt", sep="\t")  # gene name

# summary statistics
get_sumstat <- function(pheno, mode, binary) {
  setwd("~/Documents/Harpak/X inactivation/XWAS")
  if (binary) {bin <- ".glm.logistic.hybrid"} else {bin <- ".glm.linear"}
  if (mode == "ADD") {mod <- "."} else {mod <- "_dom."}
  file_name <- paste0("fe. male_X", mod, pheno, bin)
  df <- read.csv(file_name, sep="\t", colClasses = c("NULL", "integer", "character", rep("NULL",4), "character", rep("NULL",4), "numeric"))
  df <- df[!is.na(df$P),]
  df <- df[df$TEST == mode,]
}

######## XWAS ########
# variables #### USER INPUTS ####
pheno <- "K80" 
title <- "Cholelithiasis"
binary <- TRUE

df_add <- get_sumstat(pheno, "ADD", binary)
df_dom <- get_sumstat(pheno, "DOMDEV", binary)

# merge dominance and additive df with gene names and ASE data
df_xwas <- rbind(df_dom, df_add)
df_xwas <- merge(df_xwas, gene_list, by.x = "ID", by.y = "rsid") # merge with gene names
df_xwas <- merge(df_xwas, df_ase, by.x = "gene_name", by.y = "Gene") # merge with ASE
df_xwas$neg_logP <- -log(df_xwas$P, base=10)

# take the snp with highest neg log P to represent the gene
df_xwas <- df_xwas %>%
  arrange(desc(neg_logP)) %>%
  mutate(checkdup = paste0(TEST,gene_name)) 
df_xwas <- df_xwas[!duplicated(c(df_xwas$checkdup)),]

# correlation for non-zero ASE genes
df_xwas <- df_xwas %>%
  mutate(cor_test = ifelse(ASE != 0, "yes", "no")) %>%  # separate non-zero ASE
  filter(gene_name != "XIST") %>% # remove XIST as outlier
  mutate(facet_label = ifelse(TEST == "ADD", "Additive", "Dominance"))


######## PLOT ########
# Strength of association against Xi expression
ggplot(df_xwas, aes(x=ASE, y=neg_logP, color=cor_test)) +
  geom_point() +
  stat_cor(method="pearson", label.y=6, label.x=0.05, size=3.2) +
  geom_smooth(method="lm", se=F) +
  theme_classic() +
  labs(y="Strength of XWAS Association", x="Proportion of Gene Expression from Xi", title=title) +
  labs(y="", x="", title=title) + # SF
  theme(plot.title = element_text(size=13), axis.title = element_text(size=11),
       axis.text = element_text(size=8), legend.position = "none", strip.text = element_text(size = 11)) +
  theme(plot.title = element_text(size=11), axis.title = element_text(size=9), # SF
        axis.text = element_text(size=7), legend.position = "none", 
        strip.text = element_blank(), strip.background = element_blank()) +
  scale_color_manual(values=c("grey40","steelblue")) +
  geom_text_repel(aes(label=gene_name), size=3) +
  facet_wrap(vars(facet_label))


######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## TISSUE-SPECIFIC ########

# all SNPs

df_tis <- rbind(df_dom,df_add)
df_tis <- merge(df_tis, gene_list, by.x = "ID", by.y = "rsid") # merge with gene names
df_tis <- merge(df_tis, df_ase, by.x = "gene_name", by.y = "Gene") # merge with ASE
df_tis$neg_logP <- -log(df_tis$P, base=10)

# take the snp with highest neg log P to represent the gene
df_tis <- df_tis[order(df_tis$neg_logP, decreasing = TRUE),]
df_tis$checkdup <- paste0(df_tis$TEST,df_tis$gene_name,df_tis$Tissue)
df_tis <- df_tis[!duplicated(c(df_tis$checkdup)),]

# correlation for non-zero ASE genes
df_tis$cor_test <- ifelse(df_tis$ASE != 0, "yes", "no") # separate non-zero ASE
df_tis <- df_tis[df_tis$gene_name != "XIST",] # remove XIST as outlier

# make df of correlation coefficients
# cor_df <- df2 %>%
#   filter(cor_test == "yes") %>%
#   select(gene_name,TEST,Tissue,ASE,neg_logP) %>%
#   group_by(TEST, Tissue) %>%
#   summarise(correlation=cor.test(neg_logP, ASE)$estimate,
#             cor_pvalue=cor.test(neg_logP, ASE)$p.value) %>%
#   as.data.frame()
# write.table(cor_df, paste0(pheno, "_ase_bytissue.txt"), sep="\t", row.names = F, quote=F)

# rename for easier plotting
df_tis$Tissue <- replace(df_tis$Tissue, df_tis$Tissue == "Cells - EBV - transformed lymphocytes", "LCL")
df_tis$Tissue <- replace(df_tis$Tissue, df_tis$Tissue == "Skin - Not Sun Exposed (Suprapubic)", "Skin - Suprapubic")

## plot all tissues
ggplot(df_tis, aes(x=ASE, y=neg_logP, color=cor_test)) +
  geom_point(size=0.8) +
  stat_cor(method="pearson", label.y=6, label.x=0.05, size=2.5) +
  geom_smooth(method="lm", se=F, size=0.6) +
  theme_classic() +
  labs(y="-log10(P)", x="Xi to Total Expression", title=title) +
  scale_color_manual(values=c("grey40","steelblue")) +
  # geom_text_repel(aes(label=gene_name), size=1.6) +
  facet_grid(Tissue~TEST) +
  theme(legend.position = "none", plot.title = element_text(size=11), 
        axis.title = element_text(size=9), axis.text = element_text(size=7),
        strip.text.y.right = element_text(angle = 0, hjust=0), strip.background.y = element_blank(),
        strip.text = element_text(size=8))
dev.off()

## plot specific tissue
tissue <- "Stomach"   #### USER INPUT ####
df_tis_spec <- df_tis[df_tis$Tissue == tissue,]
df_tis_spec$facet_label <- ifelse(df_tis_spec$TEST == "ADD", "Additive", "Dominance")

ggplot(df_tis_spec, aes(x=ASE, y=neg_logP, color=cor_test)) +
  geom_point() +
  stat_cor(method="pearson", label.y=6, label.x=0.15) +
  geom_smooth(method="lm", se=F) +
  theme_classic() +
  labs(y="Strength of XWAS Association", x="Proportion of Gene Expression from Xi", 
       title=paste0(title, " in ", tissue)) +
  theme(plot.title = element_text(size=13), axis.title = element_text(size=11),
        axis.text = element_text(size=8), legend.position = "none", strip.text = element_text(size = 11)) +
  scale_color_manual(values=c("grey40","steelblue")) +
  geom_text_repel(aes(label=gene_name), size=3.2) +
  facet_wrap(vars(facet_label))



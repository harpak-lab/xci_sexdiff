library(ggplot2) 
library(ggrepel)
library(reshape2)
library(stringr)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyr)

setwd("~/Documents/Harpak/X inactivation/Readings/Tukiainen_nature24265-s3/")

# load ASE data for each tissue
ase_df <- read.csv("Suppl.Table.5.ASE.txt", sep="\t")
# format tissue names
ase_df <- ase_df %>% 
  mutate(Tissue = str_replace(Tissue, "ESPMCS", "Esophagus - Mucosa")) %>% mutate(Tissue = str_replace(Tissue, "PNCREAS", "Pancreas")) %>%
  mutate(Tissue = str_replace(Tissue, "STMACH", "Stomach")) %>% mutate(Tissue = str_replace(Tissue, "KDNCTX", "Kidney - Cortex")) %>%
  mutate(Tissue = str_replace(Tissue, "ESPMSL", "Esophagus - Muscularis")) %>% mutate(Tissue = str_replace(Tissue, "VAGINA", "Vagina")) %>%
  mutate(Tissue = str_replace(Tissue, "ARTAORT", "Artery - Aorta")) %>% mutate(Tissue = str_replace(Tissue, "CLNTRN", "Colon - Transverse")) %>%
  mutate(Tissue = str_replace(Tissue, "LUNG", "Lung")) %>% mutate(Tissue = str_replace(Tissue, "WHLBLD", "Whole Blood")) %>%
  mutate(Tissue = str_replace(Tissue, "ARTCRN", "Artery - Coronary")) %>% mutate(Tissue = str_replace(Tissue, "LIVER", "Liver")) %>%
  mutate(Tissue = str_replace(Tissue, "SKINS", "Skin - Not Sun Exposed (Suprapubic)")) %>% mutate(Tissue = str_replace(Tissue, "THYROID", "Thyroid")) %>%
  mutate(Tissue = str_replace(Tissue, "UTERUS", "Uterus")) %>% mutate(Tissue = str_replace(Tissue, "LCL", "Cells - EBV-transformed lymphocytes")) 
ase_df$Gene_ID <- gsub("\\..*","",ase_df$Gene_ID)
ase_df <- ase_df[c(1,2,4,10)]
write.table(ase_df, "ASE_tissue.txt", sep="\t", row.names = F, quote = F)

# load gene expression for NPX and PAR1 genes
setwd("~/Documents/Harpak/X inactivation/")
npx_df <- read.csv("meddiff_SMTSD_allgene.txt", sep="\t")
par_df <- read.csv("meddiff_SMTSD_PAR1gene.txt", sep="\t")

# y homologue
y_df <- read.csv("allNPX.txt", sep="\t")

# function to merge df
merge_df <- function(npx_df, par_df, ase_df, tis_name="sum") {
  # merge
  df <- rbind(npx_df, par_df)
  ifelse(tis_name == "sum", col_lab = c("Gene", "median"), col_lab = c("SMTSD", "Gene", "median"))
  df <- merge(df, ase_df, by="Gene") 
  # add regions
  df$Region <- ifelse(df$Gene %in% par_df$Gene, "PAR1", "NPX") 
  # add Y-homologue
  y_df <- merge(df, y_df[c(1,2)], by="Gene", all.x = T)
  y_df <- y_df %>% filter(NPY.Pair)
  # remove XIST as outlier
  df <- df[!df$Gene %in% c("XIST"),]
  return(df)
}

# ACROSS TISSUES
# get median across tissues
med_tissue <- function(df){
  colnames(df)[which(names(df) == "variable")] <- "Gene"
  df <- df %>% 
    group_by(Gene) %>%
    summarise(across(where(is.numeric), ~median(.x, na.rm=TRUE))) %>%
    as.data.frame()
  return(df)
}
ase_df_sum <- med_tissue(ase_df)
npx_df_sum <- med_tissue(npx_df)
par_df_sum <- med_tissue(par_df)

# merge df
df <- merge_df(npx_df_sum, par_df_sum, ase_df_sum)

# plot aggregate plot
ggplot(df_sum, aes(x=ASE, y=median, color=Region)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x, se=FALSE) +
  scale_color_manual(values=c("#D64B27","grey40")) +
  geom_point(data=y_df, aes(x=ASE, y=median), color="#7f171c", shape=17, size=2) +
  labs(title = "All Tissues", y="Median Female-Male TPM Difference",
       x="Proportion of Gene Expression from Xi") +
  theme_classic() +
  stat_cor(label.y=c(0.50, 0.45)) +
  theme(legend.position = "none",
        plot.title = element_text(size=13),
        axis.title = element_text(size=11),
        axis.text = element_text(size=8)) +
  geom_text_repel(aes(label=Gene), size=3.2)


# TISSUE SPECIFIC
tissue <- "Skin - Not Sun Exposed (Suprapubic)"    ## USER INPUT
ase_df_tis <- ase_df[ase_df$Tissue == tissue, c("Gene", "ASE")]
npx_df_tis <- npx_df[npx_df$SMTSD == tissue,]
par_df_tis <- par_df[par_df$SMTSD == tissue,]

# merge df
df_tissue <- merge_df(npx_df_tis, par_df_tis, ase_df_tis, tissue)

# plot tissue-specific
  ggplot(df_tissue, aes(x=ASE, y=median, color=Region)) +
    geom_point() +
    geom_smooth(method="lm", formula=y~x, se=FALSE) +
    stat_cor(label.y=c(0.67, 0.62), label.x=c(0.0,0.0), size=3.2) +
    labs(y="Median Female-Male TPM Difference",
         x="Proportion of Gene Expression from Xi", title=tissue) +
    labs(y="", x="", title = tissue) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(size=13),
          axis.title = element_text(size=11),
          axis.text = element_text(size=8)) +
    scale_color_manual(values=c("#D64B27","grey40")) +
    geom_point(data=y_df, aes(x=ASE, y=median), color="#7f171c", shape=17, size=2) +
    geom_text_repel(aes(label=Gene), size=3.2)









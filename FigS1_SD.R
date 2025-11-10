library(ggplot2) 
library(ggrepel)
library(reshape2)
library(stringr)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyr)

setwd("/Users/carriezhu/Documents/Harpak/X inactivation/XCI_code/xci_sexdiff/")
# read gene expression and ase data for NPX region
npx_df <- read.csv("meddiff_SMTSD_NPXgene.txt", sep="\t")
npx_df <- npx_df[npx_df$data == "med_diff_adj", c(1,3,4)]
ase_df <- read.csv("ASE_tissue.txt", sep="\t")

# merge dataframes
df <- merge(npx_df, ase_df, by.x=c("variable","SMTSD"), by.y=c("Gene","Tissue"))
df <- df %>% 
  dplyr::select(variable, SMTSD, value, mean_ASE)
colnames(df) <- c("Gene","SMTSD","TPMdiff", "ASE")

# get sd
sd_df <- df %>%
  filter(ASE>0) %>%
  group_by(Gene) %>%
  mutate(n=n()) %>%
  group_by(Gene, n) %>%
  summarise_at(vars(TPMdiff, ASE), sd) %>%
  as.data.frame() %>%
  na.omit() %>%
  filter(!Gene %in% c("XIST")) %>%
  arrange((TPMdiff)) 
head(sd_df)

#remove outliers
sd_df <- sd_df[sd_df$Gene != "PTCHD1",]

# plot 4x6
ggplot(sd_df, aes(x=ASE, y=TPMdiff)) +
  geom_point(color="#D64B27") +
  geom_smooth(method="lm", formula=y~x, se=FALSE, color = "#D64B27") +
  stat_cor(label.y=c(0.3), label.x=c(0.1), color="#D64B27") +
  labs(#title="Variable Escape Predicts Variation in \nSex Differences in Expression Across Tissues", 
    title = "Standard Deviation Across Tissues",
    y = "Female-Male \nGene Expression Differencess", 
    x="Proportion of Xi Expression") +
  theme_classic() +
  theme(plot.title = element_text(size=14), axis.title = element_text(size=14),
        axis.text = element_text(size=10)) +
  geom_text_repel(aes(label=Gene), size=3.2, color="#D64B27")


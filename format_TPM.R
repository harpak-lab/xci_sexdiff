library(reshape2)
library(stringr)
library(dplyr)

setwd("~/Documents/Harpak/X inactivation")

# tissue data
tissue <- read.csv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t", 
                   colClasses=c("character", rep("NULL", 5), "character", rep("NULL", 56)))
# remove brain and CML
tissue <- tissue[!str_detect(tissue$SMTSD, "Brain"),]
tissue <- tissue[!tissue$SMTSD == "Cells - Leukemia cell line (CML)",]
t_types <- unique(tissue$SMTSD)

# sex data
# 1=Male; 2=Female
sex <- read.csv("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep="\t",
                colClasses=c("character", "integer", "NULL", "NULL"))

# format TPM data to include sex and tissue metadata
format_tpm <- function(df_tpm) {
  # transpose df to have samples as rows and genes as columns
  df_tpm <- df_tpm[,c(-1)]
  df <- setNames(data.frame(t(df_tpm[,-1])), df_tpm[,1])
  df$SAMPID <- row.names(df)
  row.names(df) <- NULL
  df$SAMPID <- str_replace_all(df$SAMPID, "\\.", "-")
  # merge tissue and sex data
  df <- merge(df, tissue, by="SAMPID")
  df$SAMPID <- sub("(.*?-.*?)-.*", "\\1", df$SAMPID)
  df <- merge(df, sex, by.x = "SAMPID", by.y="SUBJID")
  return(df)
}

# NPX
df_npx <- read.csv("geneNPX_tpm_allsamples.txt", sep="\t",
              colClasses=c("character", "character", rep("numeric",17382)))
df_npx <- format_tpm(df_npx)
write.table(df, "geneNPX_TPM_filtered.txt", sep="\t", quote=FALSE)

# PAR1
df_par1 <- read.csv("genePAR1_tpm_allsamples.txt", sep="\t",
                   colClasses=c("character", "character", rep("numeric",17382)))
df_par1 <- format_tpm(df_par1)
write.table(df, "genePAR1_TPM_filtered.txt", sep="\t", quote=FALSE)





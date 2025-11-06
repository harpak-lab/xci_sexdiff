library(reshape2)
library(stringr)
library(dplyr)

med_diff <- function(df_TPM) {
  # male median to standardize
  male_median <- df_TPM %>%
    select(-c(SAMPID)) %>%
    filter(SEX == 1) %>%
    group_by(SMTSD) %>%
    summarise(across(where(is.numeric), median)) %>%
    as.data.frame()
  
  # get median by tissue and sex
  agg_df <- df_TPM %>%
    select(-c(SAMPID)) %>%
    group_by(SMTSD,SEX) %>%
    summarise(across(where(is.numeric), median)) %>%
    as.data.frame()
  
  # remove male or female only tissues
  single <- agg_df %>%
    group_by(SMTSD) %>%
    filter(n()==1) %>%
    select(SMTSD) %>%
    as.data.frame
  agg_df <- agg_df[!agg_df$SMTSD %in% single$SMTSD,]
  male_median <- male_median[!male_median$SMTSD %in% single$SMTSD,]
  
  # subtract female median - male median
  m_df <- agg_df[agg_df$SEX == 1, c(3:ncol(agg_df))]
  f_df <- agg_df[agg_df$SEX == 2, c(3:ncol(agg_df))]
  s_df <- f_df - m_df
  male_median <- male_median[,c(2:(ncol(male_median)-1))]
  # standardize by male median per tissue, per gene
  ss_adj <- s_df/male_median
  
  # add/edit labels
  agg_df <- agg_df %>%
    filter(SEX == 1) %>%
    select(SMTSD)
  ss_adj$data <- "med_diff_adj"; m_df$data <- "male_med"; f_df$data <- "female_med"
  # combine columns
  total_df <- cbind(agg_df, rbind(ss_adj, m_df, f_df))
  # make long version of dataframe
  long_df <- melt(total_df, id.vars = c("SMTSD","data"))
  long_df[sapply(long_df, is.infinite)] <- NA
  long_df <- long_df[complete.cases(long_df),]    # remove Na and Inf
  
  return(long_df)
}

df_npx <- read.csv("geneAll_TPM_filtered.txt", sep="\t")
df_npx <- med_diff(df_npx)
write.table(df_npx, "meddiff_SMTSD_NPXgene.txt", sep="\t", row.names=FALSE, quote=FALSE)

df_par1 <- read.csv("genePAR1_TPM_filtered.txt", sep="\t")
df_par1 <- med_diff(df_par1)
write.table(df_par1, "meddiff_SMTSD_PAR1gene.txt", sep="\t", row.names=FALSE, quote=FALSE)



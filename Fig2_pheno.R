setwd("/scratch/09059/xliaoyi/harpak_lab/ukb_data/fids")
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)


g43 <-  read.csv("fid41270_G43.csv", sep = ',')
# g43_report <- read.csv("fid131052.csv", sep = ',')
# g43_report <- g43_report %>% filter(X131052.0.0 != "")
# g43_flt <- g43 %>% filter(eid %in% g43_report$eid)
# g43_flt <- g43

hba1c <- read.csv("fid30750.csv", sep = ',')
hba1c_flt <- hba1c %>% select(c('eid', 'X30750.0.0')) %>% filter(!is.na(X30750.0.0)) %>% rename(hba1c = X30750.0.0)

# # remove outliers
# upper = mean(hba1c_flt$hba1c) + sd(hba1c_flt$hba1c)
# lower = mean(hba1c_flt$hba1c) - sd(hba1c_flt$hba1c)


g43_hba1c <- merge(g43, hba1c_flt, by = 'eid', all = TRUE)


# ---------------------------
# WB eids genotype file
setwd("~/Documents/Harpak/X inactivation/XCI_code")
both_gt <- read.csv("20250705_17snps_carrie_wbeids.raw", sep = '\t')

male_gt <- both_gt %>% filter(SEX == 1)
female_gt <- both_gt %>% filter(SEX == 2)

#----------------------------

# male_gt <- read.csv("../../ukb/PLINK_RECODE_OUTPUT/20250703_17snps_m.raw", sep = '\t')
# female_gt <- read.csv("../../ukb/PLINK_RECODE_OUTPUT/20250703_17snps_f.raw", sep = '\t')

# remove male code == 1
male_gt <- male_gt %>% 
  filter(rs143064101_G != 1) %>% 
  filter(rs151294283_G != 1) %>% 
  filter(rs190534439_A != 1) %>% 
  filter(rs183890070_G != 1) %>% 
  filter(rs191410333_T != 1) %>% 
  filter(rs142164313_T != 1) %>% 
  filter(rs55880880_C != 1) %>%
  filter(rs759743603_A != 1) %>% 
  filter(rs182401083_T != 1) %>% 
  filter(rs140850203_G != 1) %>% 
  filter(rs144965452_C != 1) %>% 
  filter(rs181655018_G != 1) %>% 
  filter(rs56277954_C != 1) %>% 
  filter(rs781086854_A != 1) %>% 
  filter(rs111552619_G != 1) %>% 
  filter(rs2743900_T != 1) %>% 
  filter(rs182580507_C != 1)

male_gt_pheno <- merge(male_gt, g43_hba1c, by.x = 'FID', by.y = 'eid')
female_gt_pheno <- merge(female_gt, g43_hba1c, by.x = 'FID', by.y = 'eid')


# snplist = c('rs143064101_G', 'rs151294283_G', 'rs190534439_A', 'rs183890070_G', 'rs191410333_T', 'rs142164313_T', 'rs55880880_C')
snplist = colnames(male_gt)[7:length(colnames(male_gt))]

minor_allele_list = c(
  'rs143064101_G' = "T", 
  'rs151294283_G' = "A", 
  'rs190534439_A' = "G", 
  'rs183890070_G' = "A", 
  'rs191410333_T' = "C", 
  'rs142164313_T' = "C", 
  'rs55880880_C' = "T", 
  'rs759743603_A' = "G", 
  'rs182401083_T' = "C", 
  'rs140850203_G' = "C", 
  'rs144965452_C' = "T", 
  'rs181655018_G' = "T", 
  'rs56277954_C' = "T", 
  'rs781086854_A' = "C", 
  'rs111552619_G' = "A", 
  'rs2743900_T' = "A",
  'rs182580507_C' = "T"
)
phenos = c("G43", "hba1c")


j = 1
plot <- list()
for (snp in snplist) {
  print(snp)
  minor_allele = minor_allele_list[[snp]]
  major_allele = strsplit(snp, "_")[[1]][-1]
  rsid = strsplit(snp, "_")[[1]][1]
  
  i = 1
  p <- list()
  for (pheno in phenos) {
    
    print(pheno)
    
    # process male data
    m = male_gt_pheno %>% 
      select(all_of(c("IID", snp, pheno))) %>% 
      filter(if_all(everything(), ~ !is.na(.))) %>% 
      mutate(GT = case_when(
        .data[[snp]] == 2 ~ major_allele,
        .data[[snp]] == 0 ~ minor_allele,
        TRUE      ~ NA_character_
      ) 
      ) %>% 
      mutate(GT = factor(GT, levels = c(major_allele, minor_allele)))
    
    # process female data
    f = female_gt_pheno %>% 
      select(all_of(c("IID", snp, pheno))) %>% 
      filter(if_all(everything(), ~ !is.na(.))) %>% 
      mutate(GT = case_when(
        .data[[snp]] == 2 ~ paste0(major_allele, "/", major_allele),
        .data[[snp]] == 1 ~ paste0(major_allele, "/", minor_allele),
        .data[[snp]] == 0 ~ paste0(minor_allele, "/", minor_allele),
        TRUE      ~ NA_character_
      )
      ) %>% 
      mutate(GT = factor(GT, levels = c(
        paste0(major_allele, "/", major_allele),
        paste0(major_allele, "/", minor_allele),
        paste0(minor_allele, "/", minor_allele)
      )))
    
    # ------------------------------------------------------
    # By default, REF alleles are now counted in PLINK2
    # for this SNP rs2743900_T, T is the reference allele but has lower allele frequency
    if (snp == 'rs2743900_T') {
      m = male_gt_pheno %>% 
        select(all_of(c("IID", snp, pheno))) %>% 
        filter(if_all(everything(), ~ !is.na(.))) %>% 
        mutate(GT = case_when(
          .data[[snp]] == 2 ~ major_allele,
          .data[[snp]] == 0 ~ minor_allele,
          TRUE      ~ NA_character_
        ) 
        ) %>% 
        mutate(GT = factor(GT, levels = c(minor_allele, major_allele)))
      
      # process female data
      f = female_gt_pheno %>% 
        select(all_of(c("IID", snp, pheno))) %>% 
        filter(if_all(everything(), ~ !is.na(.))) %>% 
        mutate(GT = case_when(
          .data[[snp]] == 2 ~ paste0(major_allele, "/", major_allele),
          .data[[snp]] == 1 ~ paste0(minor_allele, "/", major_allele),
          .data[[snp]] == 0 ~ paste0(minor_allele, "/", minor_allele),
          TRUE      ~ NA_character_
        )
        ) %>% 
        mutate(GT = factor(GT, levels = c(
          paste0(minor_allele, "/", minor_allele),
          paste0(minor_allele, "/", major_allele),
          paste0(major_allele, "/", major_allele)
        )))
    }
    # ------------------------------------------------------
    
    # plot
    if (pheno == 'G43') {
      # --- summary ---
      df_summary_m <- m %>%
        group_by(GT) %>%
        summarise(
          prop = mean(G43 == 1),
          n = n(),
          se = sqrt(prop * (1 - prop) / n)
        )
      
      df_summary_f <- f %>%
        group_by(GT) %>%
        summarise(
          prop = mean(G43 == 1),
          n = n(),
          se = sqrt(prop * (1 - prop) / n)
        ) 
      
      # --- flexiable same ylim ---
      all_proportion_data <- bind_rows(
        df_summary_m %>% select(prop, se),
        df_summary_f %>% select(prop, se)
      )
      
      min_val <- min(all_proportion_data$prop - all_proportion_data$se)
      max_val <- max(all_proportion_data$prop + all_proportion_data$se)
      common_ylim <- c(min_val, max_val)
      
      # --- plot ---
      pm <- ggplot(df_summary_m, aes(x = GT, y = prop)) +
        geom_point(size = 3, color = "steelblue") +
        # Add error bars based on SE
        geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "steelblue") +
        ylab("") +
        xlab(rsid) +
        ylim(common_ylim) +
        theme_bw() + 
        ggtitle("Males") +
        scale_x_discrete(labels = function(x) paste0(x, "\nn=", df_summary_m$n[match(x, df_summary_m$GT)]))
      
      pm2 <- ggplot(df_summary_m, aes(x = GT, y = prop)) +
        geom_point(size = 3, color = "steelblue") +
        # Add error bars based on SE
        geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "steelblue") +
        ylab("") +
        xlab(rsid) +
        theme_bw() + 
        ggtitle("Males Zoomed in") +
        scale_x_discrete(labels = function(x) paste0(x, "\nn=", df_summary_m$n[match(x, df_summary_m$GT)]))
      
      pf <- ggplot(df_summary_f, aes(x = GT, y = prop)) +
        geom_point(size = 3, color = "indianred") +
        # Add error bars based on SE
        geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "indianred") +
        ylab("Proportion with G43 = 1") +
        xlab(rsid) +
        ylim(common_ylim) +
        theme_bw() +
        ggtitle("Females") +
        scale_x_discrete(labels = function(x) paste0(x, "\nn=", df_summary_f$n[match(x, df_summary_f$GT)]))
      
      a <- plot_grid(pf, pm, pm2, labels = "", rel_widths = c(3, 2, 2), ncol = 3)
    } else {
      # --- summary ---
      df_summary_m <- m %>%
        group_by(GT) %>%
        summarise(
          n = n(),
          mean_hba1c = mean(hba1c),
          sd_hba1c = sd(hba1c),
          se_hba1c = sd_hba1c / sqrt(n)
        )
      
      df_summary_f <- f %>%
        group_by(GT) %>%
        summarise(
          n = n(),
          mean_hba1c = mean(hba1c),
          sd_hba1c = sd(hba1c),
          se_hba1c = sd_hba1c / sqrt(n)
        ) 
      
      # --- flexiable same ylim ---
      all_hba1c_data <- bind_rows(
        df_summary_m %>% select(mean_hba1c, se_hba1c),
        df_summary_f %>% select(mean_hba1c, se_hba1c)
      )
      
      min_val <- min(all_hba1c_data$mean_hba1c - all_hba1c_data$se_hba1c)
      max_val <- max(all_hba1c_data$mean_hba1c + all_hba1c_data$se_hba1c)
      common_ylim <- c(min_val, max_val)
      
      # --- plot ---
      
      # For Males (pm)
      pm <- ggplot(df_summary_m, aes(x = GT, y = mean_hba1c)) +
        geom_point(size = 3, color = "steelblue") + # Points for the mean
        geom_errorbar(aes(ymin = mean_hba1c-se_hba1c, ymax = mean_hba1c+se_hba1c), width = 0.2, color = "steelblue") + # Error bars for 95% CI
        ylab("") +
        xlab(rsid) +
        theme_bw() +
        ylim(common_ylim) +
        ggtitle("Males") +
        scale_x_discrete(labels = function(x) paste0(x, "\nn=", df_summary_m$n[match(x, df_summary_m$GT)]))
      
      pm2 <- ggplot(df_summary_m, aes(x = GT, y = mean_hba1c)) +
        geom_point(size = 3, color = "steelblue") + # Points for the mean
        geom_errorbar(aes(ymin = mean_hba1c-se_hba1c, ymax = mean_hba1c+se_hba1c), width = 0.2, color = "steelblue") + # Error bars for 95% CI
        ylab("") +
        xlab(rsid) +
        theme_bw() +
        ggtitle("Males Zoomed in") +
        scale_x_discrete(labels = function(x) paste0(x, "\nn=", df_summary_m$n[match(x, df_summary_m$GT)]))
      
      # For Females (pf)
      pf <- ggplot(df_summary_f, aes(x = GT, y = mean_hba1c)) +
        geom_point(size = 3, color = "indianred") + # Points for the mean
        geom_errorbar(aes(ymin = mean_hba1c-se_hba1c, ymax = mean_hba1c+se_hba1c), width = 0.2, color = "indianred") + # Error bars for 95% CI
        ylab("HbA1c") +
        xlab(rsid) +
        theme_bw() +
        ylim(common_ylim) +
        ggtitle("Females") +
        scale_x_discrete(labels = function(x) paste0(x, "\nn=", df_summary_f$n[match(x, df_summary_f$GT)]))
      
      a <- plot_grid(pf, pm, pm2, labels = "", rel_widths = c(3, 2, 2), ncol = 3)
    }
    p[[i]] <- a
    i = i + 1
  }
  plot[[j]] <- plot_grid(p[[1]], p[[2]], labels = "", rel_widths = c(1, 1))
  j = j + 1
}

allp <- plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]], plot[[6]], plot[[7]], 
                  plot[[8]], plot[[9]], plot[[10]], plot[[11]], plot[[12]], plot[[13]], plot[[14]], 
                  plot[[15]], plot[[16]], plot[[17]], 
                  labels = c(''), label_size = 12, ncol = 1)
ggsave("../../tmp/17_snps.pdf", plot = allp, width = 16, height = 37, units = "in")
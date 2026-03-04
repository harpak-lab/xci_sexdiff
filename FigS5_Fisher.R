library(reshape2)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)



# compare individual read counts using Fisher test
ase_df <- read.csv("Gylemo_ase_raw.txt", sep="\t", 
                   colClasses = c("NULL", "NULL", rep("character",4), rep("NULL", 5), rep("numeric", 4), rep("NULL",3)))

# Filter only the two individuals of interest
pair_df <- ase_df %>%
  filter(individual %in% c("UPIC", "nmXCI-1"))  # USER INPUT: [UPIC, nmXCI-1, nmXCI-2]

# create minCount vs majCount columns
pair_df <- pair_df %>%
  mutate(minCount = pmin(refCount, altCount)) %>%
  mutate(majCount = pmax(refCount, altCount))

# Function to compute Fisher test for one gene-tissue pair
fisher_gene_tissue <- function(sub) {
  if(nrow(sub) != 2) return(NA)  # confirm exactly two individuals
  mat <- matrix(c(sub$minCount, sub$majCount), nrow = 2, byrow = FALSE)
  res <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA)
  return(res)
}

# Apply across all gene × tissue pairs
results <- pair_df %>%
  group_by(gene, tissue_id, PAR) %>%
  summarise(p_value = fisher_gene_tissue(pick(everything())), .groups = "drop")

# Apply FDR correction (Bonf.)
results <- results %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

# Cap FDR at 0.05 for visualization (values above 0.05 shown as nonsignificant) and remove na
results <- results %>%
  mutate(logFDR = -log10(FDR),
         signif = ifelse(FDR < 0.05, "Significant", "Not significant")) %>%
  na.omit()

# Inspect significant results
head(results[order(results$p_value), ])
# proportion significant
nrow(results[results$signif == "Significant",])
nrow(results)
unique(results[results$signif == "Significant", c("gene")])

# max -logFDR
max_FDR <- max(results$logFDR)

# split by npx and par regions
npx_results <- results[results$PAR == "nonPAR",]
par_results <- results[results$PAR == "PAR",]

# Order genes and tissues by mean significance (for cleaner plotting)
order_results <- function(results_df) {
  gene_order <- results_df %>%
    group_by(gene) %>%
    summarise(mean_logFDR = mean(logFDR, na.rm = TRUE)) %>%
    arrange(desc(mean_logFDR)) %>%
    pull(gene)
  
  tissue_order <- results_df %>%
    group_by(tissue_id) %>%
    summarise(mean_logFDR = mean(logFDR, na.rm = TRUE)) %>%
    arrange(desc(mean_logFDR)) %>%
    pull(tissue_id)
  
  results_df$gene <- factor(results_df$gene, levels = gene_order)
  results_df$tissue_id <- factor(results_df$tissue_id, levels = tissue_order)
  
  return(results_df)
}
npx_results <- order_results(npx_results)
par_results <- order_results(par_results)


# Plot heatmap
npx <- ggplot(npx_results, aes(x = tissue_id, y = gene, fill = logFDR)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient(low = "white", high = "#D64B27", limits = c(0, max_FDR),
                      name = "-log10(FDR)",
                      na.value = "white") +
  labs(title = "Fisher's Exact Test Results by Gene and Tissue", subtitle = "nmXCI-1 and nmXCI-2") +
  theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=5),
        panel.grid = element_blank())

par <- ggplot(par_results, aes(x = tissue_id, y = gene, fill = logFDR)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient(low = "white", high = "grey10", limits = c(0, max_FDR),
                      name = "-log10(FDR)",
                      na.value = "white") +
  theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=5),
        panel.grid = element_blank())

gA <- ggplotGrob(npx) ; gB <- ggplotGrob(par)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5]) 
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth) 

plot <- grid.arrange(gA, gB, ncol=1, nrow=2, heights=c(5,1))



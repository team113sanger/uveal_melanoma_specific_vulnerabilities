library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(ggrepel)

# Load data
anno <- read.table("Data/annotation.tsv", sep = " ",
                   header = TRUE) %>% 
  select(DepMap_ID, primary_disease)

avana <- read.table("Data/avana.tsv", sep = " ",
                    header = TRUE)
colnames(avana) <- gsub("\\.", "-", colnames(avana))
avana <- na.omit(avana)

uvm_beta_scores <- read.delim("Data/MAGeCK_gene_corrected_beta.tsv")

pan_essential_genes <- read.csv("Data/pan_genes.csv")

# Tally primary diseases
primary_disease_counts <- table(anno$primary_disease)
primary_disease_counts

# Divide depmap data into SKCM and Pan-cancer
df_skin <- anno |>
  filter(primary_disease == 'Skin Cancer')

df_Nskin <- anno |>
  filter(primary_disease != 'Skin Cancer')
         
avana_sk <- avana[, intersect(df_skin$DepMap_ID, colnames(avana))]
avana_sk <- cbind(genes = avana$genes, avana_sk)

avana_nsk <- avana[, intersect(df_Nskin$DepMap_ID, colnames(avana))]
avana_nsk <- cbind(genes = avana$genes, avana_nsk)
         
cat('SKCM:', ncol(avana_sk), ';', 'Pan-Cancer:', ncol(avana_nsk), '\n')

# Select only common genes in our data and depmap
common_genes <- intersect(uvm_beta_scores$genes, avana$genes)

avana_sk_scores <- avana_sk %>% 
  filter(genes %in% common_genes)
avana_nsk_scores <- avana_nsk %>% 
  filter(genes %in% common_genes)
uvm_scores <- uvm_beta_scores %>% 
  filter(genes %in% common_genes)

# Get ranks
#uvm_scores <- uvm_beta_scores[order(uvm_beta_scores$genes),]

# need to add genes col name back...
ranked_uvm <- cbind(uvm_scores[, "genes"], apply(uvm_scores[,-1], 2, rank)) %>%
  as.data.frame() %>% 
  rename(genes = V1) %>% 
  mutate(across(-genes, as.numeric))

ranked_avana_sk <-cbind(avana_sk_scores[, "genes"], apply(avana_sk_scores[,-1], 2, rank)) %>%
  as.data.frame() %>% 
  rename(genes = V1) %>% 
  mutate(across(-genes, as.numeric))

ranked_avana_nsk <- cbind(avana_nsk_scores[, "genes"], apply(avana_nsk_scores[,-1], 2, rank)) %>%
  as.data.frame() %>% 
  rename(genes = V1) %>% 
  mutate(across(-genes, as.numeric))

# Calculate gene essentialities (mean)
uvm_ess <- uvm_scores %>%
  mutate(UVM = rowMeans(across(-genes))) %>%
  select(genes, UVM)

skcm_ess <- avana_sk_scores %>%
  mutate(SKCM = rowMeans(across(-genes))) %>%
  select(genes, SKCM)

pan_cancer_ess <- avana_nsk_scores %>%
  mutate(Pan_cancer = rowMeans(across(-genes))) %>%
  select(genes, Pan_cancer)

merged_ess <- uvm_ess %>%
  inner_join(skcm_ess, by = "genes") %>%
  inner_join(pan_cancer_ess, by = "genes")

# Calculate fold change (mann whitney test)
# test does ranking for you... do i need the previous rank step?
perform_mann_whitney <- function(data1, data2) {
  test_result <- wilcox.test(as.numeric(data1), as.numeric(data2))
  return(c(W_statistic = as.numeric(test_result$statistic), p_value = test_result$p.value))
}

# UVM vs SKCM
uvm_vs_skcm <- t(apply(cbind(ranked_uvm[, -1], ranked_avana_sk[, -1]), 1, function(row) {
  perform_mann_whitney(row[1:ncol(ranked_uvm) - 1], row[(ncol(ranked_uvm)):length(row)])
}))

uvm_vs_skcm <- as.data.frame(uvm_vs_skcm)
uvm_vs_skcm$genes <- ranked_uvm[, 1]

uvm_vs_skcm$Padj <- p.adjust(uvm_vs_skcm$p_value, method = "BH")

uvm_vs_skcm <- uvm_vs_skcm %>%
  select(genes, W_statistic, p_value, Padj)

# UVM vs Pan Cancer 
uvm_vs_pan_cancer <- t(apply(cbind(ranked_uvm[, -1], ranked_avana_nsk[, -1]), 1, function(row) {
  perform_mann_whitney(row[1:ncol(ranked_uvm) - 1], row[(ncol(ranked_uvm)):length(row)])
}))

uvm_vs_pan_cancer <- as.data.frame(uvm_vs_pan_cancer)
uvm_vs_pan_cancer$genes <- ranked_uvm[, 1]

uvm_vs_pan_cancer$Padj <- p.adjust(uvm_vs_pan_cancer$p_value, method = "BH")

uvm_vs_pan_cancer <- uvm_vs_pan_cancer %>%
  select(genes, W_statistic, p_value, Padj)

# 'LFC'...? 
uvm_medians <- apply(ranked_uvm[, -1], 1, median)
sk_medians <- apply(ranked_avana_sk[, -1], 1, median)
nsk_medians <- apply(ranked_avana_nsk[, -1], 1, median)

uvm_vs_skcm$uvm_median_rank <- uvm_medians
uvm_vs_skcm$sk_median_rank <- sk_medians
uvm_vs_skcm$FC <- uvm_vs_skcm$sk_median_rank/uvm_vs_skcm$uvm_median_rank
uvm_vs_skcm$LFC <- log2(uvm_vs_skcm$FC)
  
uvm_vs_pan_cancer$uvm_median_rank <- uvm_medians
uvm_vs_pan_cancer$nsk_median_rank <- nsk_medians
uvm_vs_pan_cancer$FC <- uvm_vs_pan_cancer$nsk_median_rank/uvm_vs_pan_cancer$uvm_median_rank
uvm_vs_pan_cancer$LFC <- log2(uvm_vs_skcm$FC)

# Combine data 
uvm_vs_skcm <- merge(uvm_vs_skcm, merged_ess[,c("genes", "UVM", "SKCM")],
                     by = "genes")
uvm_vs_pan_cancer <- merge(uvm_vs_pan_cancer, merged_ess[,c("genes", "UVM", "Pan_cancer")],
                           by = "genes")

# Pan ess genes
Log2FC_cutoff <- 1.5
padj_cutoff <- 0.05

uvm_vs_skcm_filtered <- uvm_vs_skcm |>
  filter(!genes %in% pan_essential_genes[["Gene"]]) 

uvm_vs_skcm_filtered_sig <- uvm_vs_skcm_filtered |>
  filter(LFC > Log2FC_cutoff) |>
  filter(Padj < padj_cutoff) |>
  arrange(desc(LFC))

uvm_vs_pan_cancer_filtered <- uvm_vs_pan_cancer |>
  filter(!genes %in% pan_essential_genes[["Gene"]]) 

uvm_vs_pan_cancer_filtered_sig <- uvm_vs_pan_cancer_filtered |>
  filter(LFC > Log2FC_cutoff) |>
  filter(Padj < padj_cutoff) |>
  arrange(desc(LFC))

# Box plots
top_genes <- head(uvm_vs_skcm_filtered_sig$genes, 10)

ranked_uvm_filtered <- ranked_uvm %>%
  filter(genes %in% top_genes)

ranked_avana_sk_filtered <- ranked_avana_sk %>%
  filter(genes %in% top_genes)

ranked_uvm_long <- ranked_uvm_filtered %>% 
  pivot_longer(cols = -genes, names_to = "cell_line", values_to = "Rank") %>%
  mutate(type = "UVM")

ranked_avana_sk_long <- ranked_avana_sk_filtered %>% 
  pivot_longer(cols = -genes, names_to = "cell_line", values_to = "Rank") %>%
  mutate(type = "SKCM")

box_plot_data <- bind_rows(ranked_uvm_long, ranked_avana_sk_long)
box_plot_data[["genes"]] <-
  factor(box_plot_data[["genes"]], levels = top_genes)
box_plot_data[["type"]] <-
  factor(box_plot_data[["type"]], levels = c("UVM", "SKCM"))

ggplot(box_plot_data, aes(x = genes, y = Rank, fill = type)) +
  geom_boxplot(notch = FALSE) +
 # scale_y_continuous(trans = 'log10') +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
#  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.justification = c("right", "top"),
        legend.position = c(.95, 1.1),
        panel.grid.major.y = element_line()) +
  labs(title = "Top 10 UVM-specific Gene Vulnerabilities",
       y = "Rank") +
  scale_fill_manual(values = c("UVM" = "deepskyblue4", "SKCM" = "sandybrown"))

# with stats
ggplot(box_plot_data, aes(x = type, y = Rank, fill = type)) +
  geom_boxplot(notch = FALSE) +
  geom_signif(
    comparisons = list(c("UVM", "SKCM")),
    map_signif_level = TRUE,
    test = "wilcox.test"
  ) +
  facet_wrap(~ genes, scales = "free") +
  coord_cartesian(ylim = c(0, max(box_plot_data$Rank) * 1.15)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line()) +
  labs(title = "Top 10 UVM-specific Gene Vulnerabilities",
       y = "Rank") +
  scale_fill_manual(values = c("UVM" = "deepskyblue4", "SKCM" = "sandybrown"))

# Volcano Plots 
volcano_plot_data <- uvm_vs_skcm_filtered %>%
  mutate(significant = case_when(
      Padj < 0.05 & abs(LFC) > 2 ~ ifelse(LFC > 0, "Upregulated", "Downregulated"),
      TRUE ~ "Not Significant"
    ))

top_genes <- volcano_plot_data %>% 
  filter(significant != "Not Significant") %>% 
  top_n(10, abs(LFC))

ggplot(volcano_plot_data, aes(x = LFC, y = -log10(Padj), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "indianred3", "Downregulated" = "royalblue3"),
                     guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = top_genes, aes(label = genes), size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
  labs(x = "Log2(Fold Change)", y = "-Log10(Padj)") +
  theme_bw()



  
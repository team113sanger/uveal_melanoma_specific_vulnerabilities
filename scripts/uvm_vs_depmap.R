source("src/calculate_fold_change.R")
source("src/boxplots.R")
source("src/volcano_plot.R")

library(readr)

#### Load Data ####
# Load depmap data
avana_sk_scores <- read_tsv("processed_data/avana_sk_scores.tsv")
avana_Nsk_scores <- read_tsv("processed_data/avana_Nsk_scores.tsv")

# Load UVM data
uvm_scores <- read_tsv("processed_data/uvm_scores.tsv")

# Load other
pan_essential_genes <- read_csv("data/pan_genes.csv")

#### Get ranked data ####
uvm_ranks <- rank_scores(uvm_scores)
write_tsv(uvm_ranks, "processed_data/uvm_ranks.tsv")

avana_sk_ranks <- rank_scores(avana_sk_scores)
write_tsv(avana_sk_ranks, "processed_data/avana_sk_ranks.tsv")

avana_Nsk_ranks <- rank_scores(avana_Nsk_scores)
write_tsv(avana_Nsk_ranks, "processed_data/avana_Nsk_ranks.tsv")

#### Calculate fold changes ####
# UVM vs SKCM
uvm_vs_skcm <- calculate_fold_change(
  get_mann_whitney_results(uvm_ranks, avana_sk_ranks),
  uvm_ranks,
  avana_sk_ranks
)

# UVM vs Pan Cancer
uvm_vs_pan_cancer <- calculate_fold_change(
  get_mann_whitney_results(uvm_ranks, avana_Nsk_ranks),
  uvm_ranks,
  avana_Nsk_ranks
)

# Filter out pan essential genes
uvm_vs_skcm_filtered <- filter_fc_results(
  uvm_vs_skcm, pan_essential_genes,
  signifcant_only = FALSE
)
write_tsv(uvm_vs_skcm_filtered, "processed_data/uvm_vs_skcm.tsv")

uvm_vs_pan_cancer_filtered <- filter_fc_results(
  uvm_vs_pan_cancer, pan_essential_genes,
  signifcant_only = FALSE
)
write_tsv(uvm_vs_pan_cancer_filtered, "processed_data/uvm_vs_pan_cancer.tsv")

# Filter out pan essentials and select only significant UVM specific genes
uvm_vs_skcm_filtered_sig <- filter_fc_results(
  uvm_vs_skcm, pan_essential_genes,
  signifcant_only = TRUE
)
write_tsv(uvm_vs_skcm_filtered_sig, "processed_data/uvm_vs_skcm_sig.tsv")

uvm_vs_pan_cancer_filtered_sig <- filter_fc_results(
  uvm_vs_pan_cancer, pan_essential_genes,
  signifcant_only = TRUE
)
write_tsv(uvm_vs_pan_cancer_filtered_sig, "processed_data/uvm_vs_pan_cancer_sig.tsv")

#### Box plots ####
# Get top 10 UVM specific genes
top_genes_sk <- head(uvm_vs_skcm_filtered_sig$genes, 10)
top_genes_Nsk <- head(uvm_vs_pan_cancer_filtered_sig$genes, 10)

# Create plots
uvm_vs_skcm_p <- plot_boxplot(
  prepare_boxplot_data(uvm_ranks, top_genes_sk, "UVM"),
  prepare_boxplot_data(avana_sk_ranks, top_genes_sk, "SKCM"),
  "UVM", "SKCM"
)
ggsave("plots/uvm_vs_skcm_boxplot.pdf", uvm_vs_skcm_p,
  width = 9, height = 5
)

uvm_vs_pan_p <- plot_boxplot(
  prepare_boxplot_data(uvm_ranks, top_genes_Nsk, "UVM"),
  prepare_boxplot_data(avana_Nsk_ranks, top_genes_Nsk, "Pan_cancer"),
  "UVM", "Pan_cancer"
)
ggsave("plots/uvm_vs_pan_cancer_boxplot.pdf", uvm_vs_pan_p,
  width = 9, height = 5
)

# Create plots with stats bars
uvm_vs_skcm_p <- plot_stats_boxplots(
  prepare_boxplot_data(uvm_ranks, top_genes_sk, "UVM"),
  prepare_boxplot_data(avana_sk_ranks, top_genes_sk, "SKCM"),
  "UVM", "SKCM"
)
ggsave("plots/uvm_vs_skcm_stats_boxplot.pdf", uvm_vs_skcm_p,
  width = 10, height = 8
)

uvm_vs_pan_p <- plot_stats_boxplots(
  prepare_boxplot_data(uvm_ranks, top_genes_Nsk, "UVM"),
  prepare_boxplot_data(avana_Nsk_ranks, top_genes_Nsk, "Pan_cancer"),
  "UVM", "Pan_cancer"
)
ggsave("plots/uvm_vs_pan_cancer_stats_boxplot.pdf", uvm_vs_pan_p,
  width = 10, height = 8
)

#### Volcano Plots ####
# Create plots
uvm_vs_skcm_p <- plot_volcano(prepare_volcano_data(uvm_vs_skcm_filtered))
ggsave("plots/uvm_vs_skcm_volcano.pdf", uvm_vs_skcm_p,
  width = 9, height = 5
)

uvm_vs_pan_cancer_p <- plot_volcano(prepare_volcano_data(uvm_vs_pan_cancer_filtered))
ggsave("plots/uvm_vs_pan_cancer_volcano.pdf", uvm_vs_pan_cancer_p,
  width = 9, height = 5
)

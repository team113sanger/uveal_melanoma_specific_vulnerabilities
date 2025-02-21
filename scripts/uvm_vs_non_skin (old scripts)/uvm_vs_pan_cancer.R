source("src/calculate_fold_change.R")
source("src/boxplots.R")
source("src/volcano_plot.R")

library(readr)
library(purrr)

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
avana_sk_ranks <- rank_scores(avana_sk_scores)
avana_Nsk_ranks <- rank_scores(avana_Nsk_scores)

# Save files
walk2(
  list(uvm_ranks, avana_sk_ranks, avana_Nsk_ranks),
  list(
    "processed_data/uvm_ranks.tsv",
    "processed_data/avana_sk_ranks.tsv",
    "processed_data/avana_Nsk_ranks.tsv"
  ), write_tsv
)

#### Calculate fold changes ####
# UVM vs SKCM
uvm_vs_skcm <- get_mann_whitney_results(uvm_ranks, avana_sk_ranks) |>
  calculate_fold_change(stats_df = _, uvm_ranks, avana_sk_ranks)

# UVM vs Pan Cancer
uvm_vs_pan_cancer <- get_mann_whitney_results(uvm_ranks, avana_Nsk_ranks) |>
  calculate_fold_change(stats_df = _, uvm_ranks, avana_Nsk_ranks)

# Filter out pan essential genes
filter_fc_and_write <- function(df, genes_to_remove, sig_only, file_name) {
  filtered_fcs <- filter_fc_results(
    df = df,
    genes_to_remove = genes_to_remove,
    signifcant_only = sig_only
  )
  write_tsv(filtered_fcs, file.path("processed_data", file_name))
  return(filtered_fcs)
}

filtered_fold_change_dfs <- pmap(
  list(
    df = list(
      uvm_vs_skcm,
      uvm_vs_pan_cancer,
      uvm_vs_skcm,
      uvm_vs_pan_cancer
    ),
    genes_to_remove = list(rep(pan_essential_genes, 4)),
    sig_only = list(FALSE, FALSE, TRUE, TRUE),
    file_name = list(
      "uvm_vs_skcm.tsv",
      "uvm_vs_pan_cancer.tsv",
      "uvm_vs_skcm_sig.tsv",
      "uvm_vs_pan_cancer_sig.tsv"
    )
  ),
  .f = filter_fc_and_write
)

uvm_vs_skcm_filtered <- filtered_fold_change_dfs[[1]]
uvm_vs_pan_cancer_filtered <- filtered_fold_change_dfs[[2]]
uvm_vs_skcm_filtered_sig <- filtered_fold_change_dfs[[3]]
uvm_vs_pan_cancer_filtered_sig <- filtered_fold_change_dfs[[4]]

#### Box plots ####
# Get top 10 UVM specific genes
top_genes_sk <- head(uvm_vs_skcm_filtered_sig[["genes"]], 10)
top_genes_nUVM <- head(uvm_vs_pan_cancer_filtered_sig[["genes"]], 10)

# Create plots
create_and_save_boxplots <- function(
  ranks_df_a, ranks_df_b, top_genes, label_a,
  label_b, stats = FALSE, file_name
) {
  plot_df <- prepare_boxplot_data(
    ranks_df_a, ranks_df_b, top_genes, label_a, label_b
  )
  
  plot <- plot_boxplot(plot_df = plot_df, stats = stats)
  
  ggsave(
    file.path("plots", file_name),
    plot = plot,
    width = ifelse(stats, 10, 9),
    height = ifelse(stats, 8, 5)
  )
}

params <- list(
  ranks_df_a = list(uvm_ranks, uvm_ranks, uvm_ranks, uvm_ranks),
  ranks_df_b = list(
  avana_sk_ranks,
  avana_nUVM_ranks,
  avana_sk_ranks,
  avana_nUVM_ranks
  ),
  top_genes = list(top_genes_sk, top_genes_nUVM, top_genes_sk, top_genes_nUVM),
  label_a = list("UVM", "UVM", "UVM", "UVM"),
  label_b = list("SKCM", "Pan_cancer", "SKCM", "Pan_cancer"),
  stats = list(FALSE, FALSE, TRUE, TRUE),
  file_name = list(
  "uvm_vs_skcm_boxplot.pdf",
  "uvm_vs_pan_cancer_boxplot.pdf",
  "uvm_vs_skcm_stats_boxplot.pdf",
  "uvm_vs_pan_cancer_stats_boxplot.pdf"
)
)

pmap(params, create_and_save_boxplots)

#### Volcano Plots ####
create_and_save_volcano <- function(df, file_name) {
  label_significant_genes(df) |>
    plot_volcano(df = _) |>
    ggsave(
      file.path("plots", file_name),
      plot = _,
      width = 9, height = 5
    )
}

walk2(list(uvm_vs_skcm_filtered, uvm_vs_pan_cancer_filtered),
  list("uvm_vs_skcm_volcano.pdf", "uvm_vs_pan_cancer_volcano.pdf"),
  .f = create_and_save_volcano
)

source("src/calculate_fold_change.R")
source("src/boxplots.R")
source("src/volcano_plot.R")

library(readr)
library(purrr)
library(dplyr)

#### Load Data ####
# Load depmap scores
pan_cancer_scores <- read_tsv("processed_data/pan_cancer_scores.tsv")

# Load UVM scores
uvm_scores <- read_tsv("processed_data/uvm_scores.tsv")

# Load other
pan_essential_genes <- read_csv("data/pan_genes.csv")

#### Get ranked data ####
uvm_ranks <- rank_scores(uvm_scores)
pan_cancer_ranks <- rank_scores(pan_cancer_scores)

# Save files
walk2(
  list(uvm_ranks, pan_cancer_ranks),
  list(
    "results/uvm_ranks.tsv",
    "results/pan_cancer_ranks.tsv"
  ), write_tsv
)

#### Calculate fold changes ####
# UVM vs Pan Cancer
uvm_vs_pan_cancer <- get_mann_whitney_results(uvm_ranks, pan_cancer_ranks) |>
  calculate_fold_change(stats_df = _, uvm_ranks, pan_cancer_ranks)

# Filter out pan essential genes
filter_fc_and_write <- function(df, genes_to_remove, sig_only, file_name) {
  filtered_fcs <- filter_fc_results(
    df = df,
    genes_to_remove = genes_to_remove,
    signifcant_only = sig_only
  )
  write_tsv(filtered_fcs, file.path("results", file_name))
  return(filtered_fcs)
}

filtered_fold_change_dfs <- pmap(
  list(
    df = list(
      uvm_vs_pan_cancer,
      uvm_vs_pan_cancer
    ),
    genes_to_remove = list(rep(pan_essential_genes, 2)),
    sig_only = list(FALSE, TRUE),
    file_name = list(
      "uvm_vs_pan_cancer.tsv",
      "uvm_vs_pan_cancer_sig.tsv"
    )
  ),
  .f = filter_fc_and_write
)

uvm_vs_pan_cancer_filtered <- filtered_fold_change_dfs[[1]]
uvm_vs_pan_cancer_filtered_sig <- filtered_fold_change_dfs[[2]]

#### Box plots ####
# Get top 10 UVM specific genes
top_genes <- head(uvm_vs_pan_cancer_filtered_sig[["genes"]], 10)
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
  ranks_df_a = rep(list(uvm_ranks), 2),
  ranks_df_b = rep(list(pan_cancer_ranks), 2),
  top_genes = rep(list(top_genes), 2),
  label_a = list("UVM", "UVM"),
  label_b = list("Pan_cancer", "Pan_cancer"),
  stats = list(FALSE, TRUE),
  file_name = list(
  "uvm_vs_pan_cancer_boxplot.pdf",
  "uvm_vs_pan_cancer_stats_boxplot.pdf"
)
)

pmap(params, create_and_save_boxplots)

#### Volcano Plots ####
create_and_save_volcano <- function(df, file_name) {
  df <- label_significant_genes(df)
  plot <- plot_volcano(df)

  ggsave(file.path("plots", paste0(file_name, ".pdf")), plot, width = 8, height = 7)
  ggsave(file.path("plots", paste0(file_name, ".png")), plot, width = 8, height = 7, dpi = 300)
}

create_and_save_volcano(uvm_vs_pan_cancer_filtered, "uvm_vs_pan_cancer_volcano")

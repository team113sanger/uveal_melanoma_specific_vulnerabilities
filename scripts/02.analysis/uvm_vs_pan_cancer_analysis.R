source("scripts/02.analysis/calculate_fold_change.R")
source("scripts/03.visualisation/boxplots.R")

library(readr)
library(purrr)
library(dplyr)
library(tidyr)

#### Load Data ####
# Load depmap scores
pan_cancer_scores <- read_tsv("data/processed/pan_cancer_scores.tsv")

# Load UVM scores
uvm_scores <- read_tsv("data/processed/uvm_scores.tsv")

# Load other
pan_essential_genes <- read_csv("data/metadata/pan_genes.csv")

#### Get ranked data ####
uvm_ranks <- rank_scores(uvm_scores)
pan_cancer_ranks <- rank_scores(pan_cancer_scores)

# Save files
walk2(
  list(uvm_ranks, pan_cancer_ranks),
  list(
    "results/tables/uvm_ranks.tsv",
    "results/tables/pan_cancer_ranks.tsv"
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
  write_tsv(filtered_fcs, file.path("results", "tables", file_name))
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

# Process ranks dfs in prep for plotting
uvm_vs_pan_cancer_sig <- filtered_fold_change_dfs[[2]]

top_genes <- head(uvm_vs_pan_cancer_sig[["genes"]], 10)

processed_ranks <- prepare_boxplot_data(
    uvm_ranks, pan_cancer_ranks, top_genes, "UVM", "Pan-Cancer"
  )
write_tsv(processed_ranks, "results/tables/top_10_uvm_pan_cancer_ranks.tsv")
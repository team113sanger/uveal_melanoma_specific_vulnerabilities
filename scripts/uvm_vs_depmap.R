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

# #### Get ranked data ####
uvm_ranks <- rank_scores(uvm_scores)

write_tsv(uvm_ranks, "processed_data/uvm_ranks.tsv")

avana_sk_ranks <- rank_scores(avana_sk_scores)
write_tsv(avana_sk_ranks, "processed_data/avana_sk_ranks.tsv")

avana_Nsk_ranks <- rank_scores(avana_Nsk_scores)
write_tsv(avana_Nsk_ranks, "processed_data/avana_Nsk_ranks.tsv")

# Could do something like this:
# map(list(uvm_scores, avana_sk_scores, avana_Nsk_scores), rank_scores) |> 
# walk2(.x =_, 
#       .y = list("processed_data/uvm_ranks.tsv", 
#                 "processed_data/avana_sk_ranks.tsv",
#                 "processed_data/avana_Nsk_ranks.tsv"), 
#        .f = ~write_tsv(.x,.y))


#### Calculate fold changes ####
# UVM vs SKCM
uvm_vs_skcm <- get_mann_whitney_results(uvm_ranks, avana_sk_ranks) |> 
                calculate_fold_change(stats_df = _, uvm_ranks,avana_sk_ranks)

# UVM vs Pan Cancer
uvm_vs_pan_cancer <- get_mann_whitney_results(uvm_ranks, avana_Nsk_ranks) |>
  calculate_fold_change(stats_df = _, uvm_ranks, avana_Nsk_ranks)

# Filter out pan essential genes
uvm_vs_skcm_filtered <- filter_fc_results(
  uvm_vs_skcm, pan_essential_genes,
  signifcant_only = FALSE
)
write_tsv(uvm_vs_skcm_filtered, "processed_data/uvm_vs_skcm.tsv")

filter_fc_and_write <- function(dataset, genes, sig_only, file_path){
    filtered_fcs <- filter_fc_results(df = dataset, 
                                      genes_to_remove = genes, 
                                      signifcant_only = sig_only)
    write_tsv(filtered_fcs, file = file_path)
    return(filtered_fcs)
}

fold_change_datasets <- pmap(list(
    dataset = list(uvm_vs_skcm, 
                   uvm_vs_pan_cancer,
                   uvm_vs_skcm, 
                   uvm_vs_pan_cancer),
    genes = list(rep(pan_essential_genes, 4)),
    sig_only = list(FALSE, FALSE, TRUE, TRUE),
    file_path = list("processed_data/uvm_vs_skcm.tsv",
                     "processed_data/uvm_vs_pan_cancer.tsv",
                     "processed_data/uvm_vs_skcm_sig.tsv",
                     "processed_data/uvm_vs_pan_cancer_sig.tsv")
                     ),
    .f = filter_fc_and_write)

uvm_vs_skcm_filtered <- fold_change_datasets[[1]]
uvm_vs_pan_cancer_filtered <- fold_change_datasets[[2]]
uvm_vs_skcm_filtered_sig <- fold_change_datasets[[3]]
uvm_vs_pan_cancer_filtered_sig <- fold_change_datasets[[4]]

#### Box plots ####
# Get top 10 UVM specific genes
top_genes_sk <- head(uvm_vs_skcm_filtered_sig[["genes"]], 10)
top_genes_Nsk <- head(uvm_vs_pan_cancer_filtered_sig[["genes"]], 10)

# Create plots
prepare_boxplot_data(
    uvm_ranks, avana_sk_ranks, top_genes_sk, "UVM", "SKCM"
  ) |>
plot_boxplot(plot_df = _) |>
ggsave(
  "plots/uvm_vs_skcm_boxplot.pdf", plot = _,
  width = 9, height = 5
)

prepare_boxplot_data(
    uvm_ranks, avana_Nsk_ranks, top_genes_Nsk, "UVM", "Pan_cancer"
  ) |>
plot_boxplot(plot_df = _
) |>
ggsave(
  "plots/uvm_vs_pan_cancer_boxplot.pdf", plot = _,
  width = 9, height = 5
)

# Create plots with stats bars
prepare_boxplot_data(
    uvm_ranks, avana_sk_ranks, top_genes_sk, "UVM", "SKCM") |>
plot_stats_boxplots(plot_df = _
) |>
ggsave(
  "plots/uvm_vs_skcm_stats_boxplot.pdf", plot = _,
  width = 10, height = 8
)

prepare_boxplot_data(
    uvm_ranks, avana_Nsk_ranks, top_genes_Nsk, "UVM", "Pan_cancer") |>
plot_stats_boxplots(plot_df = _) |>
ggsave(
  "plots/uvm_vs_pan_cancer_stats_boxplot.pdf", plot = _,
  width = 10, height = 8)




#### Volcano Plots ####

# label_significant_genes(uvm_vs_skcm_filtered) |>
# plot_volcano(df = _) |>
# ggsave(
#   "plots/uvm_vs_skcm_volcano.pdf", plot = _,
#   width = 9, height = 5
# )

# label_significant_genes(uvm_vs_pan_cancer_filtered) |>
# plot_volcano(df = _) |>
# ggsave(
#   "plots/uvm_vs_pan_cancer_volcano.pdf", plot = _,
#   width = 9, height = 5
# )


volcano_sig_genes <- function(dataset, filepath){
label_significant_genes(dataset) |>
plot_volcano(df = _) |>
ggsave(
  filepath, plot = _,
  width = 9, height = 5
)
}

walk2(list(uvm_vs_skcm_filtered,uvm_vs_pan_cancer_filtered),
    list("plots/uvm_vs_skcm_volcano.pdf", "plots/uvm_vs_pan_cancer_volcano.pdf"), 
    .f = volcano_sig_genes)


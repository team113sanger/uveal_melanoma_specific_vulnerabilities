library(dplyr)
library(readr)

# Load depmap data
avana_sk_ranks <- read_tsv("processed_data/avana_sk_ranks.tsv")
avana_Nsk_ranks <- read_tsv("processed_data/avana_Nsk_ranks.tsv")

# Load UVM data
uvm_ranks <- read_tsv("processed_data/uvm_ranks.tsv")

# Load pan essential genes list
pan_essential_genes <- read_csv("raw_data/pan_genes.csv")

# Perform wilcox test
perform_mann_whitney <- function(data_a, data_b) {
  test_result <- wilcox.test(as.numeric(data_a), as.numeric(data_b))
  return(c(
    W_statistic = as.numeric(test_result[["statistic"]]),
    p_value = test_result[["p.value"]]
  ))
}

get_mann_whitney_results <- function(ranks_df, ranks_df_2) {
  df <- t(apply(cbind(ranks_df[, -1], ranks_df_2[, -1]), 1, function(row) {
    perform_mann_whitney(
      row[1:ncol(ranks_df) - 1], row[(ncol(ranks_df)):length(row)]
    )
  }))

  df <- as.data.frame(df) %>%
    mutate(genes = ranks_df[["genes"]]) %>%
    mutate(Padj = p.adjust(p_value, method = "BH")) %>%
    select(genes, W_statistic, p_value, Padj)

  return(df)
}

# UVM vs SKCM
uvm_vs_skcm <- get_mann_whitney_results(uvm_ranks, avana_sk_ranks)

# UVM vs Pan Cancer
uvm_vs_pan_cancer <- get_mann_whitney_results(uvm_ranks, avana_Nsk_ranks)

# Add median ranks and calculate fold changes
calculate_fold_change <- function(df, ranks_df, ranks_df_2) {
  new_df <- df %>%
  mutate(
    median_rank_a = apply(ranks_df[, -1], 1, median),
    median_rank_b = apply(ranks_df_2[, -1], 1, median),
    FC = median_rank_b / median_rank_a,
    LFC = log2(FC)
  )

  return(new_df)
}

uvm_vs_skcm <- calculate_fold_change(uvm_vs_skcm, uvm_ranks, avana_sk_ranks)
uvm_vs_pan_cancer <- calculate_fold_change(uvm_vs_pan_cancer, uvm_ranks, avana_Nsk_ranks)

filter_fc_results <- function(df, genes_to_remove, signifcant_only = FALSE) {
  df_filtered <- df %>%
    filter(!genes %in% genes_to_remove[["Gene"]])

  if(signifcant_only) {
    Log2FC_cutoff <- 1.5
    padj_cutoff <- 0.05

    df_filtered_sig <- df_filtered %>%
      filter(LFC > Log2FC_cutoff & Padj < padj_cutoff) %>%
      arrange(desc(LFC))

    return(df_filtered_sig)
  } else {
    return(df_filtered)
  }
}

# Filter out pan essential genes
uvm_vs_skcm_filtered <- filter_fc_results(
  uvm_vs_skcm, pan_essential_genes, signifcant_only = FALSE
  )
write_tsv(uvm_vs_skcm_filtered, "processed_data/uvm_vs_skcm.tsv")

uvm_vs_pan_cancer_filtered <- filter_fc_results(
  uvm_vs_pan_cancer, pan_essential_genes, signifcant_only = FALSE
  )
write_tsv(uvm_vs_pan_cancer_filtered, "processed_data/uvm_vs_pan_cancer.tsv")

# Filter out pan essentials and select only significant genes
uvm_vs_skcm_filtered_sig <- filter_fc_results(
  uvm_vs_skcm_filtered, pan_essential_genes, signifcant_only = TRUE)
write_tsv(uvm_vs_skcm_filtered_sig, "processed_data/uvm_vs_skcm_sig.tsv")

uvm_vs_pan_cancer_filtered_sig <- filter_fc_results(
  uvm_vs_pan_cancer_filtered, pan_essential_genes, signifcant_only = TRUE)
write_tsv(uvm_vs_pan_cancer_filtered_sig, "processed_data/uvm_vs_pan_cancer_sig.tsv")

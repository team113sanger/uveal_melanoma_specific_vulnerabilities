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
    mutate(genes =  ranks_df[["genes"]]) %>%
    mutate(Padj = p.adjust(p_value, method = "BH")) %>%
    select(genes, W_statistic, p_value, Padj)

  return(df)
}

# UVM vs SKCM
uvm_vs_skcm <- get_mann_whitney_results(ranked_uvm, ranked_avana_sk)

# UVM vs Pan Cancer
uvm_vs_pan_cancer <- get_mann_whitney_results(ranked_uvm, ranked_avana_Nsk)

# Add median ranks and calculate fold changes
uvm_vs_skcm <- uvm_vs_skcm %>%
  mutate(
    uvm_median_rank = apply(ranked_uvm[, -1], 1, median),
    sk_median_rank = apply(ranked_avana_sk[, -1], 1, median),
    FC = sk_median_rank / uvm_median_rank,
    LFC = log2(FC)
  )

uvm_vs_pan_cancer <- uvm_vs_pan_cancer %>%
  mutate(
    uvm_median_rank = apply(ranked_uvm[, -1], 1, median),
    Nsk_median_rank = apply(ranked_avana_Nsk[, -1], 1, median),
    FC = Nsk_median_rank / uvm_median_rank,
    LFC = log2(FC)
  )

# Filter out pan ess genes
uvm_vs_skcm_filtered <- uvm_vs_skcm |>
  filter(!genes %in% pan_essential_genes[["Gene"]])
write_tsv(uvm_vs_skcm_filtered, "processed_data/uvm_vs_skcm.tsv")

uvm_vs_pan_cancer_filtered <- uvm_vs_pan_cancer |>
  filter(!genes %in% pan_essential_genes[["Gene"]])
write_tsv(uvm_vs_pan_cancer_filtered, "processed_data/uvm_vs_pan_cancer.tsv")

# Filter for significant results
Log2FC_cutoff <- 1.5
padj_cutoff <- 0.05

uvm_vs_skcm_filtered_sig <- uvm_vs_skcm_filtered |>
  filter(LFC > Log2FC_cutoff) |>
  filter(Padj < padj_cutoff) |>
  arrange(desc(LFC))
write_tsv(uvm_vs_skcm_filtered_sig, "processed_data/uvm_vs_skcm_sig.tsv")

uvm_vs_pan_cancer_filtered_sig <- uvm_vs_pan_cancer_filtered |>
  filter(LFC > Log2FC_cutoff) |>
  filter(Padj < padj_cutoff) |>
  arrange(desc(LFC))
write_tsv(uvm_vs_pan_cancer_filtered_sig, "processed_data/uvm_vs_pan_cancer_sig.tsv")

library(dplyr)

rank_scores <- function(df) {
  ranked_df <- cbind(df[, "genes"], apply(df[, -1], 2, rank)) %>%
    as.data.frame() %>%
    mutate(across(-genes, as.numeric))

  return(ranked_df)
}

perform_mann_whitney <- function(data_a, data_b) {
  test_result <- wilcox.test(as.numeric(data_a), as.numeric(data_b))
  return(c(
    W_statistic = as.numeric(test_result[["statistic"]]),
    p_value = test_result[["p.value"]]
  ))
}

# Apply mann whitney test to each row across two dataframes of ranks
get_mann_whitney_results <- function(ranks_df, ranks_df_2) {
  df <- t(apply(cbind(ranks_df[, -1], ranks_df_2[, -1]), 1, function(row) {
    perform_mann_whitney(
      row[1:ncol(ranks_df) - 1], row[(ncol(ranks_df)):length(row)]
    )
  }))

  stats_df <- as.data.frame(df) %>%
    mutate(genes = ranks_df[["genes"]]) %>%
    mutate(Padj = p.adjust(p_value, method = "BH")) %>%
    select(genes, W_statistic, p_value, Padj)

  return(stats_df)
}

calculate_fold_change <- function(stats_df, ranks_df, ranks_df_2) {
  new_df <- stats_df %>%
    mutate(
      median_rank_a = apply(ranks_df[, -1], 1, median),
      median_rank_b = apply(ranks_df_2[, -1], 1, median),
      FC = median_rank_b / median_rank_a,
      log2FC = log2(FC)
    )

  return(new_df)
}

filter_fc_results <- function(df, genes_to_remove, signifcant_only = FALSE) {
  df_filtered <- df %>%
    filter(!genes %in% genes_to_remove[["Gene"]])

  if (signifcant_only) {
    Log2FC_cutoff <- 1.80
    padj_cutoff <- 0.01

    df_filtered_sig <- df_filtered %>%
      filter(log2FC > Log2FC_cutoff & Padj < padj_cutoff) %>%
      arrange(desc(log2FC))

    return(df_filtered_sig)
  } else {
    return(df_filtered |>
      arrange(desc(log2FC)))
  }
}

calc_avg_gene_essentiality <- function(scores_df) {
  scores_df |>
    mutate(mean_essentiality = rowMeans(across(-genes), na.rm = TRUE)) |>
    select(genes, mean_essentiality)
}

add_gene_essentiality <- function(scores_df_a, scores_df_b, lfc_results_df) {
  gene_ess_a <- calc_avg_gene_essentiality(scores_df_a)
  gene_ess_b <- calc_avg_gene_essentiality(scores_df_b)

  results_df <- lfc_results_df |>
    left_join(gene_ess_a, by = "genes") |>
    rename(mean_essentiality_a = mean_essentiality) |>
    left_join(gene_ess_b, by = "genes") |>
    rename(mean_essentiality_b = mean_essentiality)

  return(results_df)
}


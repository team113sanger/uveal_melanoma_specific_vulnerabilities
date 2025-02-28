library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# Convert dataframes into long format and combine
process_ranks_df <- function(ranks_df, top_genes, label) {
  filtered_df <- ranks_df |>
    filter(genes %in% top_genes) |>
    arrange(match(genes, top_genes)) |>
    rename(gene = genes)

  df_long <- filtered_df |>
    pivot_longer(
      cols = -gene, names_to = "cell_line", values_to = "rank"
    ) |>
    mutate(group = label)

  return(df_long)
}

prepare_boxplot_data <- function(
    ranks_df_a, ranks_df_b, top_genes, label_a, label_b) {
  df_long_a <- process_ranks_df(ranks_df_a, top_genes, label_a)
  df_long_b <- process_ranks_df(ranks_df_b, top_genes, label_b)

  combined_df <- bind_rows(df_long_a, df_long_b)

  combined_df[["gene"]] <-
    factor(combined_df[["gene"]], levels = top_genes)

  combined_df[["group"]] <-
    factor(combined_df[["group"]], levels = c(label_a, label_b))

  return(combined_df)
}

plot_boxplot <- function(plot_df, stats = FALSE) {
  labels <- unique(plot_df[["group"]])

  fill_colors <- setNames(
    c("#226E9C", "#AB1866"), labels
  )

  if (stats) {
    aes_mapping <- aes(x = group, y = rank, fill = group)
  } else {
    aes_mapping <- aes(x = gene, y = rank, fill = group)
  }

  p <- ggplot(plot_df, aes_mapping) +
    geom_boxplot(
      notch = FALSE,
      outlier.shape = 21,
      outlier.color = "black",
      outlier.fill = NA
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.63, size = 12, face = "italic"),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.justification = c("right", "top"),
      panel.grid.major.y = element_line(),
      axis.ticks = element_line(color = "black"),
    ) +
    labs(title = "Gene Rank Comparison", x = "Gene", y = "Rank") +
    scale_fill_manual(values = fill_colors)

  if (stats) {
    p +
      facet_wrap(~gene, scales = "free") +
      geom_signif(
        comparisons = list(as.character(labels)),
        map_signif_level = TRUE,
        test = "wilcox.test"
      ) +
      coord_cartesian(
        ylim = c(0, max(plot_df[["rank"]]) * 1.15)
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "plain"),
      ) +
      scale_x_discrete(guide = guide_axis(angle = 0))
  } else {
    p + theme(legend.justification = c("right", "top"))
  }
}

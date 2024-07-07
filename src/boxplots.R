library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

prepare_boxplot_data <- function(ranks_df, top_genes, label) {
  filtered_df <- ranks_df %>%
    filter(genes %in% top_genes)

  df_long <- filtered_df %>%
    pivot_longer(
      cols = -genes, names_to = "cell_line", values_to = "rank"
    ) %>%
    mutate(type = label)

  df_long[["genes"]] <-
    factor(df_long[["genes"]], levels = top_genes)

  return(df_long)
}

plot_boxplot <- function(df_a, df_b, label_a, label_b) {
  plot_df <- bind_rows(df_a, df_b)

  plot_df[["type"]] <-
    factor(plot_df[["type"]], levels = c(label_a, label_b))

  fill_colors <- setNames(
    c("deepskyblue4", "sandybrown"), c(label_a, label_b)
  )

  ggplot(plot_df, aes(x = genes, y = rank, fill = type)) +
    geom_boxplot(notch = FALSE) +
    # scale_y_continuous(trans = 'log10') +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    #  scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.justification = c("right", "top"),
      legend.position.inside = c(.95, 1.1),
      panel.grid.major.y = element_line()
    ) +
    labs(
      title = paste0("Top 10 ", label_a, "-specific Gene Vulnerabilities"),
      y = "Rank"
    ) +
    scale_fill_manual(values = fill_colors)
}

plot_stats_boxplots <- function(df_a, df_b, label_a, label_b) {
  plot_df <- bind_rows(df_a, df_b)

  plot_df[["type"]] <-
    factor(plot_df[["type"]], levels = c(label_a, label_b))

  fill_colors <- setNames(
    c("deepskyblue4", "sandybrown"), c(label_a, label_b)
  )

  ggplot(plot_df, aes(x = type, y = rank, fill = type)) +
    geom_boxplot(notch = FALSE) +
    geom_signif(
      comparisons = list(c(label_a, label_b)),
      map_signif_level = TRUE,
      test = "wilcox.test"
    ) +
    facet_wrap(~genes, scales = "free") +
    coord_cartesian(ylim = c(0, max(plot_df[["rank"]]) * 1.15)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.position = "none",
      panel.grid.major.y = element_line()
    ) +
    labs(
      title = paste0("Top 10 ", label_a, "-specific Gene Vulnerabilities"),
      y = "Rank"
    ) +
    scale_fill_manual(values = fill_colors)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# Convert dataframes into long format and combine 
prepare_boxplot_data <- function(
    ranks_df_a, ranks_df_b, top_genes, label_a, label_b) {
    # This function is defined within another function
    # You could try defining it outside 
    # Would make this make sense to have it 
    # accesible outside of the parent function
  process_ranks_df <- function(ranks_df, label) { # This function is defined within another function
    filtered_df <- ranks_df |>
      filter(genes %in% top_genes)

    df_long <- filtered_df |>
      pivot_longer(
        cols = -genes, names_to = "cell_line", values_to = "rank"
      ) |>
      mutate(type = label)

    return(df_long)
  }

  df_long_a <- process_ranks_df(ranks_df_a, label_a)
  df_long_b <- process_ranks_df(ranks_df_b, label_b)

  combined_df <- bind_rows(df_long_a, df_long_b)

  combined_df[["genes"]] <-
    factor(combined_df[["genes"]], levels = top_genes)

  combined_df[["type"]] <-
    factor(combined_df[["type"]], levels = c(label_a, label_b))

  return(combined_df)
}

plot_boxplot <- function(plot_df) {
  labels <- unique(plot_df[["type"]])

  fill_colors <- setNames(
    c("deepskyblue4", "sandybrown"), labels
  )

  ggplot(plot_df, aes(x = genes, y = rank, fill = type)) +
    geom_boxplot(
      notch = FALSE,
      outlier.shape = 21,
      outlier.color = "black",
      outlier.fill = NA
    ) +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.justification = c("right", "top"),
      # legend.position = c(.95, 1.1),
      panel.grid.major.y = element_line()
    ) +
    labs(
      title = paste0("Top 10 ", labels[1], "-specific Gene Vulnerabilities"),
      y = "Rank"
    ) +
    scale_fill_manual(values = fill_colors)
}
# How much of plot boxplots and plot_stats boxplots is shared? 
# Could you think of a way to add a "plot-type" argument 
# and then add logical statement to do different things when 
# one type of plot is true? 


plot_stats_boxplots <- function(plot_df, label_a, label_b) {
  labels <- unique(plot_df[["type"]])

  fill_colors <- setNames(
    c("deepskyblue4", "sandybrown"), labels
  )

  ggplot(plot_df, aes(x = type, y = rank, fill = type)) +
    geom_boxplot(
      notch = FALSE,
      outlier.shape = 21,
      outlier.color = "black",
      outlier.fill = NA
    ) +
    geom_signif(
      comparisons = list(as.character(labels)),
      map_signif_level = TRUE,
      test = "wilcox.test"
    ) +
    facet_wrap(~genes, scales = "free") +
    coord_cartesian(
      ylim = c(0, max(plot_df[["rank"]]) * 1.15)
    ) +
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
      title = paste0("Top 10 ", labels[1], "-specific Gene Vulnerabilities"),
      y = "Rank"
    ) +
    scale_fill_manual(values = fill_colors)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(ggrepel)

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

plot_boxplot <- function(data_a, data_b, label_a, label_b) {
  plot_data <- bind_rows(data_a, data_b)

  plot_data[["type"]] <-
    factor(plot_data[["type"]], levels = c(label_a, label_b))

  fill_colors <- setNames(
    c("deepskyblue4", "sandybrown"), c(label_a, label_b)
  )

  ggplot(plot_data, aes(x = genes, y = rank, fill = type)) +
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
      title = "Top 10 UVM-specific Gene Vulnerabilities",
      y = "Rank"
    ) +
    scale_fill_manual(values = fill_colors)
}

plot_stats_boxplots <- function(data_a, data_b, label_a, label_b) {
  plot_data <- bind_rows(data_a, data_b)

  plot_data[["type"]] <-
    factor(plot_data[["type"]], levels = c(label_a, label_b))

  fill_colors <- setNames(
    c("deepskyblue4", "sandybrown"), c(label_a, label_b)
  )

  ggplot(plot_data, aes(x = type, y = rank, fill = type)) +
    geom_boxplot(notch = FALSE) +
    geom_signif(
      comparisons = list(c(label_a, label_b)),
      map_signif_level = TRUE,
      test = "wilcox.test"
    ) +
    facet_wrap(~genes, scales = "free") +
    coord_cartesian(ylim = c(0, max(plot_data[["rank"]]) * 1.15)) +
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
      title = "Top 10 UVM-specific Gene Vulnerabilities",
      y = "Rank"
    ) +
    scale_fill_manual(values = fill_colors)
}

# Label significant genes
prepare_volcano_data <- function(df) {
  plot_data <- df %>%
    mutate(significant = case_when(
      Padj < 0.05 & abs(LFC) > 2 ~ ifelse(LFC > 0, "Upregulated", "Downregulated"),
      TRUE ~ "Not Significant"
    ))
}

plot_volcano <- function(df) {
  top_genes <- df %>%
    filter(significant != "Not Significant") %>%
    top_n(10, abs(LFC))

  ggplot(df, aes(x = LFC, y = -log10(Padj), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Not Significant" = "gray", "Upregulated" = "indianred3", "Downregulated" = "royalblue3"
      ),
    guide = "none"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    data = top_genes, aes(label = genes), size = 3,
           box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
  labs(x = "Log2(Fold Change)", y = "-Log10(Padj)") +
  theme_bw()
}
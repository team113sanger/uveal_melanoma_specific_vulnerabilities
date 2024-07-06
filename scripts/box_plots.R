library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# Load data
uvm_ranks <- read_tsv("processed_data/uvm_ranks.tsv")
avana_sk_ranks <- read_tsv("processed_data/avana_sk_ranks.tsv")
avana_Nsk_ranks <- read_tsv("processed_data/avana_Nsk_ranks.tsv")
uvm_vs_skcm_sig <- read_tsv("processed_data/uvm_vs_skcm_sig.tsv")
uvm_vs_pan_cancer_sig <- read_tsv("processed_data/uvm_vs_pan_cancer_sig.tsv")

# Get top genes
top_genes_sk <- head(uvm_vs_skcm_sig$genes, 10)
top_genes_Nsk <- head(uvm_vs_pan_cancer_sig$genes, 10)

#
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

uvm_ranks_plot_data <- prepare_boxplot_data(uvm_ranks, top_genes_sk, "UVM")
avana_sk_ranks_plot_data <- prepare_boxplot_data(avana_sk_ranks, top_genes_sk, "SKCM")

uvm_ranks_plot_data_2 <- prepare_boxplot_data(uvm_ranks, top_genes_Nsk, "UVM")
avana_Nsk_ranks_plot_data <- prepare_boxplot_data(avana_Nsk_ranks, top_genes_Nsk, "Pan_cancer")

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

uvm_vs_skcm_p <- plot_boxplot(uvm_ranks_plot_data, avana_sk_ranks_plot_data,
                              "UVM", "SKCM")
ggsave("plots/uvm_vs_skcm_boxplot.pdf", uvm_vs_skcm_p,
  width = 8, height = 5
)

uvm_vs_pan_p <- plot_boxplot(uvm_ranks_plot_data_2, avana_Nsk_ranks_plot_data,
                              "UVM", "Pan_cancer")
ggsave("plots/uvm_vs_pan_cancer_boxplot.pdf", uvm_vs_pan_p,
  width = 8, height = 5
)

# with stats
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

uvm_vs_skcm_p <- plot_stats_boxplots(uvm_ranks_plot_data, avana_sk_ranks_plot_data,
                                    "UVM", "SKCM")
ggsave("plots/uvm_vs_skcm_stats_boxplot.pdf", uvm_vs_skcm_p,
  width = 10, height = 7
)

uvm_vs_pan_p <- plot_stats_boxplots(uvm_ranks_plot_data_2, avana_Nsk_ranks_plot_data,
                                    "UVM", "Pan_cancer")
ggsave("plots/uvm_vs_pan_cancer_stats_boxplot.pdf", uvm_vs_pan_p,
  width = 10, height = 7
)
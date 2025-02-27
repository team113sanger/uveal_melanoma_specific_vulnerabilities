source("scripts/03.visualisation/boxplots.R")
source("scripts/03.visualisation/volcano_plot.R")

library(readr)
library(purrr)
library(dplyr)

uvm_vs_pan_cancer <- read_tsv("results/tables/uvm_vs_pan_cancer.tsv")
uvm_vs_pan_cancer_sig <- read_tsv("results/tables/uvm_vs_pan_cancer_sig.tsv")

#### Box plots ####
# Get top 10 UVM specific genes
top_genes <- head(uvm_vs_pan_cancer_sig[["genes"]], 10)
# Create plots
create_and_save_boxplots <- function(
  ranks_df_a, ranks_df_b, top_genes, label_a,
  label_b, stats = FALSE, file_name
) {
  plot_df <- prepare_boxplot_data(
    ranks_df_a, ranks_df_b, top_genes, label_a, label_b
  )
  
  plot <- plot_boxplot(plot_df = plot_df, stats = stats)
  
  ggsave(
    file.path("results", "figures", file_name),
    plot = plot,
    width = ifelse(stats, 10, 9),
    height = ifelse(stats, 8, 5)
  )
}

params <- list(
  ranks_df_a = rep(list(uvm_ranks), 2),
  ranks_df_b = rep(list(pan_cancer_ranks), 2),
  top_genes = rep(list(top_genes), 2),
  label_a = list("UVM", "UVM"),
  label_b = list("Pan_cancer", "Pan_cancer"),
  stats = list(FALSE, TRUE),
  file_name = list(
  "uvm_vs_pan_cancer_boxplot.pdf",
  "uvm_vs_pan_cancer_stats_boxplot.pdf"
)
)

pmap(params, create_and_save_boxplots)

#### Volcano Plots ####
create_and_save_volcano <- function(df, file_name) {
  df <- label_significant_genes(df)
  plot <- plot_volcano(df)

  ggsave(file.path("results", "figures", paste0(file_name, ".pdf")), plot, width = 8, height = 7)
  ggsave(file.path("results", "figures", paste0(file_name, ".png")), plot, width = 8, height = 7, dpi = 300)
}

create_and_save_volcano(uvm_vs_pan_cancer_filtered, "uvm_vs_pan_cancer_volcano")

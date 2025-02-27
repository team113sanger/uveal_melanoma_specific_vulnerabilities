source("scripts/03.visualisation/boxplots.R")
source("scripts/03.visualisation/volcano_plot.R")

library(readr)
library(purrr)
library(dplyr)

top_ranks <- read_tsv("results/tables/top_10_uvm_pan_cancer_ranks.tsv")
uvm_vs_pan_cancer <- read_tsv("results/tables/uvm_vs_pan_cancer.tsv")
uvm_vs_pan_cancer_sig <- read_tsv("results/tables/uvm_vs_pan_cancer_sig.tsv")

# Add levels
top_ranks[["gene"]] <-
    factor(top_ranks[["gene"]], levels = unique(top_ranks$gene))

top_ranks[["group"]] <-
    factor(top_ranks[["group"]], levels = c("UVM", "Pan-Cancer"))

#### Box plots ####
# Create plots
create_and_save_boxplots <- function(plot_df, stats, file_name){
    plot <- plot_boxplot(plot_df = plot_df, stats = stats)
  
    ggsave(
        file.path("results", "figures", paste0(file_name, ".pdf")),
        plot = plot,
        width = ifelse(stats, 10, 9),
        height = ifelse(stats, 9, 7)
    )
    ggsave(
        file.path("results", "figures", paste0(file_name, ".png")),
        plot = plot,
        width = ifelse(stats, 10, 9),
        height = ifelse(stats, 9, 7),
        dpi = 300
    )
}

params <- list(
  plot_df = rep(list(top_ranks), 2),
  stats = list(FALSE, TRUE),
  file_name = list(
  "uvm_vs_pan_cancer_boxplot",
  "uvm_vs_pan_cancer_stats_boxplot"
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

create_and_save_volcano(uvm_vs_pan_cancer, "uvm_vs_pan_cancer_volcano")

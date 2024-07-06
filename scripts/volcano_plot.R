library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load data
uvm_vs_skcm <- read_tsv("processed_data/uvm_vs_skcm.tsv")
uvm_vs_pan_cancer <- read_tsv("processed_data/uvm_vs_pan_cancer.tsv")

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

uvm_vs_skcm_p <- plot_volcano(prepare_volcano_data(uvm_vs_skcm))
ggsave("plots/uvm_vs_skcm_volcano.pdf", uvm_vs_skcm_p,
       width = 9, height = 5)

uvm_vs_pan_cancer_p <- plot_volcano(prepare_volcano_data(uvm_vs_pan_cancer))
ggsave("plots/uvm_vs_pan_cancer_volcano.pdf", uvm_vs_pan_cancer_p,
       width = 9, height = 5)
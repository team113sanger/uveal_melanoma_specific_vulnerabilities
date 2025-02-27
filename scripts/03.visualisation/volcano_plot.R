library(dplyr)
library(ggplot2)
library(ggrepel)

label_significant_genes <- function(df) {
  df %>%
    mutate(
      significant = case_when(
        Padj < 0.01 & abs(log2FC) > 1.8 ~ ifelse(log2FC > 0, "Upregulated", "Downregulated"),
        TRUE ~ "Not Significant"
      )
    )
}
plot_volcano <- function(df) {
  top_genes <- df %>%
    filter(significant != "Not Significant") %>%
    arrange(desc(log2FC)) %>%
    head(10)
    
  ggplot(df, aes(x = log2FC, y = -log10(Padj), color = significant)) +
      geom_point(shape = 19, alpha = 0.85) +
      scale_color_manual(
        values = c("Not Significant" = "#808080", "Upregulated" = "#226E9C", "Downregulated" = "#2E8B57")
      ) +
      geom_vline(xintercept = c(-1.8, 1.8), linetype = "dashed", color = "black", alpha = 0.5) +
      # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_text_repel(
        data = top_genes, aes(label = genes),
        color = "black", fontface = "italic", size = 4,
        min.segment.length = 0
      ) +
      labs(
        x = expression('log'[2]~'(Fold change)'),
        y = expression('log'[10]~'(1 / Padj)')
      ) +
      theme_classic() +
      theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 11, color = "black"),
        legend.position = "none"  
      ) +
      scale_x_reverse() 
}

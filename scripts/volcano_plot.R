# Volcano Plots
volcano_plot_data <- uvm_vs_skcm_filtered %>%
  mutate(significant = case_when(
    Padj < 0.05 & abs(LFC) > 2 ~ ifelse(LFC > 0, "Upregulated", "Downregulated"),
    TRUE ~ "Not Significant"
  ))

top_genes <- volcano_plot_data %>%
  filter(significant != "Not Significant") %>%
  top_n(10, abs(LFC))

ggplot(volcano_plot_data, aes(x = LFC, y = -log10(Padj), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = c("Not Significant" = "gray", "Upregulated" = "indianred3", "Downregulated" = "royalblue3"),
    guide = "none"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = top_genes, aes(label = genes), size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
  labs(x = "Log2(Fold Change)", y = "-Log10(Padj)") +
  theme_bw()

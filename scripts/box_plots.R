# Box plots
top_genes <- head(uvm_vs_skcm_filtered_sig$genes, 10)

ranked_uvm_filtered <- ranked_uvm %>%
  filter(genes %in% top_genes)

ranked_avana_sk_filtered <- ranked_avana_sk %>%
  filter(genes %in% top_genes)

ranked_uvm_long <- ranked_uvm_filtered %>%
  pivot_longer(cols = -genes, names_to = "cell_line", values_to = "Rank") %>%
  mutate(type = "UVM")

ranked_avana_sk_long <- ranked_avana_sk_filtered %>%
  pivot_longer(cols = -genes, names_to = "cell_line", values_to = "Rank") %>%
  mutate(type = "SKCM")

box_plot_data <- bind_rows(ranked_uvm_long, ranked_avana_sk_long)
box_plot_data[["genes"]] <-
  factor(box_plot_data[["genes"]], levels = top_genes)
box_plot_data[["type"]] <-
  factor(box_plot_data[["type"]], levels = c("UVM", "SKCM"))

ggplot(box_plot_data, aes(x = genes, y = Rank, fill = type)) +
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
    legend.position = c(.95, 1.1),
    panel.grid.major.y = element_line()
  ) +
  labs(
    title = "Top 10 UVM-specific Gene Vulnerabilities",
    y = "Rank"
  ) +
  scale_fill_manual(values = c("UVM" = "deepskyblue4", "SKCM" = "sandybrown"))

# with stats
ggplot(box_plot_data, aes(x = type, y = Rank, fill = type)) +
  geom_boxplot(notch = FALSE) +
  geom_signif(
    comparisons = list(c("UVM", "SKCM")),
    map_signif_level = TRUE,
    test = "wilcox.test"
  ) +
  facet_wrap(~genes, scales = "free") +
  coord_cartesian(ylim = c(0, max(box_plot_data$Rank) * 1.15)) +
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
  scale_fill_manual(values = c("UVM" = "deepskyblue4", "SKCM" = "sandybrown"))

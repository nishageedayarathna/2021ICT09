
# Stage 6: Volcano Plot
# Alopecia Areata Gene Expression Analysis

# Load required package
library(ggplot2)

# Load differential expression results
de_results <- read.csv(
  "results/tables/differential_expression_results.csv"
)

# Create significance category
de_results$Significance <- ifelse(
  de_results$adj.P.Val < 0.05,
  "Significant",
  "Not Significant"
)

# Create volcano plot
volcano_plot <- ggplot(
  de_results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, 1e-300)),
    color = Significance
  )
) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot of Differential Gene Expression",
    subtitle = "AAP Lesional Samples vs Normal Controls",
    x = "Log Fold Change (logFC)",
    y = "-Log10 Adjusted P-value",
    color = "Expression Result"
  ) +
  theme_minimal()

# Display plot
print(volcano_plot)

# Save plot
ggsave(
  "results/figures/volcano_plot.png",
  plot = volcano_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# PROJECT: Alopecia Areata Bioinformatics Analysis
# STAGE 2: Data Cleaning and Quality Control
# DATASET: GSE68801
# ============================================================

# 1. Set project directory
setwd("D:/Alopecia_Areata_Bioinformatics")

# 2. Load required packages
library(Biobase)
library(ggplot2)
library(dplyr)

# 3. Load saved expression matrix and metadata
expression_data <- readRDS(
  "data/raw/GSE68801_expression.rds"
)

sample_info <- read.csv(
  "data/raw/GSE68801_sample_metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# 4. Display dimensions
cat("\nExpression matrix dimensions:\n")
print(dim(expression_data))

# 5. Check missing values
missing_values <- sum(is.na(expression_data))

cat("\nTotal missing expression values:\n")
print(missing_values)

# 6. Check infinite values
infinite_values <- sum(is.infinite(expression_data))

cat("\nTotal infinite values:\n")
print(infinite_values)

# 7. Check duplicate probe IDs
duplicate_probes <- sum(
  duplicated(rownames(expression_data))
)

cat("\nDuplicate probe IDs:\n")
print(duplicate_probes)

# 8. Check expression value range
cat("\nExpression value summary:\n")
print(summary(as.vector(expression_data)))

# 9. Check sample metadata dimensions
cat("\nSample metadata dimensions:\n")
print(dim(sample_info))

# 10. Check missing metadata values
cat("\nMissing metadata values:\n")
print(colSums(is.na(sample_info)))

# 11. Calculate expression values per sample
sample_means <- colMeans(
  expression_data,
  na.rm = TRUE
)

sample_medians <- apply(
  expression_data,
  2,
  median,
  na.rm = TRUE
)

sample_sd <- apply(
  expression_data,
  2,
  sd,
  na.rm = TRUE
)

# 12. Create sample quality summary
sample_qc <- data.frame(
  Sample = colnames(expression_data),
  Mean = sample_means,
  Median = sample_medians,
  Standard_Deviation = sample_sd
)

# 13. Save sample quality results
write.csv(
  sample_qc,
  "results/tables/sample_quality_summary.csv",
  row.names = FALSE
)

# 14. Create boxplot of sample expression distributions
png(
  "results/figures/expression_boxplot.png",
  width = 1400,
  height = 900,
  res = 150
)

boxplot(
  expression_data,
  outline = FALSE,
  las = 2,
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Expression Values"
)

dev.off()

# 15. Display final quality-control summary
cat("\n========== QUALITY CONTROL SUMMARY ==========\n")

cat("Missing expression values:", missing_values, "\n")
cat("Infinite expression values:", infinite_values, "\n")
cat("Duplicate probe IDs:", duplicate_probes, "\n")

cat("\nQuality control checks completed.\n")

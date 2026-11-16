
# ============================================================
# PROJECT: Alopecia Areata Bioinformatics Analysis
# STAGE 3: Exploratory Data Analysis
# PART 1: Inspect Sample Information
# DATASET: GSE68801
# ============================================================

# 1. Set project directory
setwd("D:/Alopecia_Areata_Bioinformatics")

# 2. Load required packages
library(ggplot2)
library(dplyr)

# 3. Load expression data
expression_data <- readRDS(
  "data/raw/GSE68801_expression.rds"
)

# 4. Load sample metadata
sample_info <- read.csv(
  "data/raw/GSE68801_sample_metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# 5. Check dimensions
cat("\nExpression matrix dimensions:\n")
print(dim(expression_data))

cat("\nSample metadata dimensions:\n")
print(dim(sample_info))

# 6. Check sample alignment
cat("\nAre expression samples and metadata aligned?\n")

print(
  identical(
    colnames(expression_data),
    rownames(sample_info)
  )
)

# 7. Display sample titles
cat("\n========== SAMPLE TITLES ==========\n")

print(sample_info$title)

# 8. Identify sample characteristic columns
characteristic_columns <- grep(
  "characteristics_ch1",
  colnames(sample_info),
  value = TRUE,
  ignore.case = TRUE
)

# 9. Display all sample characteristics
cat("\n========== SAMPLE CHARACTERISTICS ==========\n")

for (column in characteristic_columns) {
  
  cat("\nColumn:", column, "\n")
  
  print(
    sort(
      table(
        as.character(sample_info[[column]]),
        useNA = "ifany"
      ),
      decreasing = TRUE
    )
  )
  
}

# 10. Save metadata inspection results
sink("results/tables/eda_sample_metadata_summary.txt")

cat("GSE68801 SAMPLE METADATA INSPECTION\n\n")

cat("Expression dimensions:\n")
print(dim(expression_data))

cat("\nMetadata dimensions:\n")
print(dim(sample_info))

cat("\nSample titles:\n")
print(sample_info$title)

cat("\nSample characteristics:\n")

for (column in characteristic_columns) {
  
  cat("\n", column, "\n")
  
  print(
    table(
      as.character(sample_info[[column]]),
      useNA = "ifany"
    )
  )
  
}

sink()

cat("\nSample metadata inspection completed!\n")


# ============================================================
# STAGE 3: DEFINE RESEARCH SAMPLE GROUPS
# ============================================================

# Extract disease status
disease_status <- sub(
  "disease status: ",
  "",
  sample_info$characteristics_ch1.1
)

# Extract skin type
skin_type <- sub(
  "skintype \\(l=lesional, nl=nonlesional, nc=normal control\\): ",
  "",
  sample_info$characteristics_ch1.2
)

# Create sample group information
sample_groups <- data.frame(
  Sample_ID = rownames(sample_info),
  Sample_Title = sample_info$title,
  Disease_Status = disease_status,
  Skin_Type = skin_type,
  stringsAsFactors = FALSE
)

# Define primary comparison groups
sample_groups$Research_Group <- ifelse(
  disease_status == "AAP" & skin_type == "L",
  "AAP_Lesional",
  ifelse(
    disease_status == "Normal" & skin_type == "NC",
    "Normal_Control",
    "Excluded"
  )
)

# Display group counts
cat("\nPrimary research group counts:\n")

print(
  table(sample_groups$Research_Group)
)

# Select the 58 samples
selected_samples <- sample_groups[
  sample_groups$Research_Group != "Excluded",
]

# Extract corresponding expression data
selected_expression <- expression_data[
  ,
  selected_samples$Sample_ID,
  drop = FALSE
]

# Save sample group information
write.csv(
  sample_groups,
  "results/tables/sample_group_information.csv",
  row.names = FALSE
)

# Save selected sample information
write.csv(
  selected_samples,
  "results/tables/selected_research_samples.csv",
  row.names = FALSE
)

# Verify selected data dimensions
cat("\nSelected expression dimensions:\n")
print(dim(selected_expression))

cat("\nSelected group counts:\n")
print(table(selected_samples$Research_Group))


# ============================================================
# EXPRESSION DISTRIBUTION VISUALIZATION
# ============================================================

# Convert expression matrix into long format
expression_df <- data.frame(
  Expression = as.vector(selected_expression),
  Sample = rep(
    colnames(selected_expression),
    each = nrow(selected_expression)
  )
)

# Add sample group information
group_lookup <- selected_samples[
  ,
  c("Sample_ID", "Research_Group")
]

expression_df$Research_Group <- group_lookup$Research_Group[
  match(expression_df$Sample, group_lookup$Sample_ID)
]

# Create distribution plot
expression_plot <- ggplot(
  expression_df,
  aes(
    x = Research_Group,
    y = Expression,
    fill = Research_Group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  labs(
    title = "Gene Expression Distribution",
    subtitle = "AAP Lesional Samples vs Normal Controls",
    x = "Sample Group",
    y = "Expression Value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 15,
      hjust = 1
    )
  )

# Display plot
print(expression_plot)

# Save plot
ggsave(
  "results/figures/group_expression_distribution.png",
  plot = expression_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ==========================================
# PCA - Principal Component Analysis
# ==========================================

# Calculate variance for each gene probe
gene_variances <- apply(selected_expression, 1, var)

# Select the 500 most variable probes
top_variable_probes <- names(
  sort(gene_variances, decreasing = TRUE)
)[1:500]

# Create PCA input (samples as rows, probes as columns)
pca_input <- t(
  selected_expression[top_variable_probes, ]
)

# Perform PCA
pca_result <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)

# Calculate percentage of variance explained
variance_explained <- (
  pca_result$sdev^2 / sum(pca_result$sdev^2)
) * 100

# Create PCA data frame
pca_data <- data.frame(
  Sample_ID = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# Add research group information
pca_data$Research_Group <- selected_samples$Research_Group[
  match(pca_data$Sample_ID, selected_samples$Sample_ID)
]

# Draw PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Research_Group)
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Gene Expression Profiles",
    subtitle = "AAP Lesional Samples vs Normal Controls",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    color = "Sample Group"
  ) +
  theme_minimal()

print(pca_plot)

# Save PCA plot
ggsave(
  "results/figures/pca_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Save PCA variance information
write.csv(
  data.frame(
    Principal_Component = paste0("PC", 1:10),
    Variance_Explained = variance_explained[1:10]
  ),
  "results/tables/pca_variance_summary.csv",
  row.names = FALSE
)

cat("PCA analysis completed successfully!\n")


# ==========================================
# Check Microarray Batch Distribution
# ==========================================

# Extract microarray batch information
batch_info <- sub(
  "^microarray batch: *",
  "",
  sample_info$characteristics_ch1.6
)

# Create metadata table
batch_summary <- data.frame(
  Sample_ID = rownames(sample_info),
  Research_Group = sample_groups$Research_Group,
  Microarray_Batch = batch_info,
  stringsAsFactors = FALSE
)

# Keep only selected research samples
selected_batch_summary <- batch_summary[
  batch_summary$Research_Group != "Excluded",
]

# Display batch distribution
batch_table <- table(
  selected_batch_summary$Research_Group,
  selected_batch_summary$Microarray_Batch
)

print(batch_table)

# Save batch information
write.csv(
  selected_batch_summary,
  "results/tables/selected_sample_batch_information.csv",
  row.names = FALSE
)

cat("Batch distribution check completed!\n")



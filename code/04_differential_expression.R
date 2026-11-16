
# ==========================================
# Stage 4: Differential Gene Expression
# ==========================================

# Set project directory
setwd("D:/Alopecia_Areata_Bioinformatics")

# Load required packages
library(limma)

# Load expression data
expression_data <- readRDS(
  "data/raw/GSE68801_expression.rds"
)

# Load sample metadata
sample_info <- read.csv(
  "data/raw/GSE68801_sample_metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Extract disease status and skin type
disease_status <- sub(
  "disease status: ",
  "",
  sample_info$characteristics_ch1.1
)

skin_type <- sub(
  "skintype \\(l=lesional, nl=nonlesional, nc=normal control\\): ",
  "",
  sample_info$characteristics_ch1.2
)

# Extract microarray batch
batch_info <- sub(
  "^microarray batch: *",
  "",
  sample_info$characteristics_ch1.6
)

# Create sample information table
sample_groups <- data.frame(
  Sample_ID = rownames(sample_info),
  Disease_Status = disease_status,
  Skin_Type = skin_type,
  Microarray_Batch = batch_info,
  stringsAsFactors = FALSE
)

# Define research groups
sample_groups$Research_Group <- ifelse(
  disease_status == "AAP" & skin_type == "L",
  "AAP_Lesional",
  ifelse(
    disease_status == "Normal" & skin_type == "NC",
    "Normal_Control",
    "Excluded"
  )
)

# Select only the two research groups
selected_samples <- sample_groups[
  sample_groups$Research_Group != "Excluded",
]

# Match expression matrix to selected samples
selected_expression <- expression_data[
  ,
  selected_samples$Sample_ID,
  drop = FALSE
]

# Set reference group
selected_samples$Research_Group <- factor(
  selected_samples$Research_Group,
  levels = c("Normal_Control", "AAP_Lesional")
)

# Convert batch to factor
selected_samples$Microarray_Batch <- factor(
  selected_samples$Microarray_Batch
)

# Check sample alignment
stopifnot(
  identical(
    colnames(selected_expression),
    selected_samples$Sample_ID
  )
)

# Display group and batch information
print(table(selected_samples$Research_Group))
print(table(
  selected_samples$Research_Group,
  selected_samples$Microarray_Batch
))

# ==========================================
# Create Statistical Model
# ==========================================

# Include both research group and microarray batch
design <- model.matrix(
  ~ Microarray_Batch + Research_Group,
  data = selected_samples
)

print(design)
print(dim(design))

# ==========================================
# Run limma Analysis
# ==========================================

# Fit linear model
fit <- lmFit(
  selected_expression,
  design
)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract differential expression results
results <- topTable(
  fit,
  coef = "Research_GroupAAP_Lesional",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Add probe IDs as a column
results$Probe_ID <- rownames(results)

# Reorder columns
results <- results[
  ,
  c(
    "Probe_ID",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]

# ==========================================
# Identify Significant Probes
# ==========================================

significant_results <- results[
  results$adj.P.Val < 0.05,
]

# Display results
cat("Total probes analysed:", nrow(results), "\n")

cat(
  "Significant probes (adjusted p-value < 0.05):",
  nrow(significant_results),
  "\n"
)

# Display top 10 results
print(head(results, 10))

# ==========================================
# Save Results
# ==========================================

write.csv(
  results,
  "results/tables/differential_expression_results.csv",
  row.names = FALSE
)

write.csv(
  significant_results,
  "results/tables/significant_differential_expression.csv",
  row.names = FALSE
)

# Save model information
write.csv(
  data.frame(
    Sample_ID = selected_samples$Sample_ID,
    Research_Group = selected_samples$Research_Group,
    Microarray_Batch = selected_samples$Microarray_Batch
  ),
  "results/tables/differential_expression_sample_info.csv",
  row.names = FALSE
)

cat("Differential expression analysis completed!\n")

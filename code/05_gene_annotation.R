
# ==========================================
# Stage 5: Gene Annotation
# Project: Alopecia Areata Gene Expression Analysis
# ==========================================

setwd("D:/Alopecia_Areata_Bioinformatics")

# Load required packages
library(AnnotationDbi)
library(hgu133plus2.db)

# ------------------------------------------
# 1. Load differential expression results
# ------------------------------------------

de_results <- read.csv(
  "results/tables/differential_expression_results.csv",
  stringsAsFactors = FALSE
)

cat("Differential expression results loaded.\n")
cat("Total probes:", nrow(de_results), "\n\n")


# ------------------------------------------
# 2. Map probe IDs to gene information
# ------------------------------------------

annotation_results <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = unique(de_results$Probe_ID),
  keytype = "PROBEID",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

cat("Probe annotation completed.\n")


# ------------------------------------------
# 3. Merge annotation with DE results
# ------------------------------------------

annotated_results <- merge(
  de_results,
  annotation_results,
  by.x = "Probe_ID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

# Save complete annotated results
write.csv(
  annotated_results,
  "results/tables/annotated_differential_expression.csv",
  row.names = FALSE
)


# ------------------------------------------
# 4. Select significant annotated probes
# ------------------------------------------

significant_annotated <- annotated_results[
  !is.na(annotated_results$adj.P.Val) &
    annotated_results$adj.P.Val < 0.05,
]

write.csv(
  significant_annotated,
  "results/tables/significant_annotated_genes.csv",
  row.names = FALSE
)


# ------------------------------------------
# 5. Calculate annotation summary
# ------------------------------------------

total_probes <- length(unique(de_results$Probe_ID))

annotated_probes <- length(
  unique(annotation_results$PROBEID[
    !is.na(annotation_results$SYMBOL)
  ])
)

significant_probes <- length(
  unique(significant_annotated$Probe_ID)
)

unique_significant_genes <- length(
  unique(na.omit(significant_annotated$SYMBOL))
)

cat("\n========== ANNOTATION SUMMARY ==========\n")

cat("Total probes analysed:", total_probes, "\n")
cat("Probes with gene symbols:", annotated_probes, "\n")
cat("Significant probes:", significant_probes, "\n")
cat("Unique significant gene symbols:",
    unique_significant_genes, "\n")


# ------------------------------------------
# 6. Display top 15 significant annotations
# ------------------------------------------

top_15_genes <- significant_annotated[
  order(significant_annotated$adj.P.Val),
]

cat("\n========== TOP 15 SIGNIFICANT RESULTS ==========\n")

print(
  head(
    top_15_genes[
      ,
      c("Probe_ID", "SYMBOL", "GENENAME",
        "logFC", "P.Value", "adj.P.Val")
    ],
    15
  )
)


# ------------------------------------------
# 7. Save annotation summary
# ------------------------------------------

summary_text <- c(
  "GENE ANNOTATION SUMMARY",
  "========================",
  paste("Total probes:", total_probes),
  paste("Probes with gene symbols:", annotated_probes),
  paste("Significant probes:", significant_probes),
  paste("Unique significant gene symbols:",
        unique_significant_genes)
)

writeLines(
  summary_text,
  "results/tables/gene_annotation_summary.txt"
)

cat("\nGene annotation completed successfully!\n")

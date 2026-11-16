
# ============================================================
# PROJECT: Alopecia Areata Bioinformatics Analysis
# STAGE 1: Download and Import GEO Dataset
# DATASET: GSE68801
# SOURCE: NCBI Gene Expression Omnibus (GEO)
# ============================================================


# ------------------------------------------------------------
# 1. SET PROJECT DIRECTORY
# ------------------------------------------------------------

setwd("D:/Alopecia_Areata_Bioinformatics")

cat("Project directory:", getwd(), "\n")


# ------------------------------------------------------------
# 2. LOAD REQUIRED PACKAGE
# ------------------------------------------------------------

library(GEOquery)
library(Biobase)

cat("Required packages loaded successfully.\n")


# ------------------------------------------------------------
# 3. CONFIGURE DOWNLOAD SETTINGS
# ------------------------------------------------------------

options(timeout = 600)

options(download.file.method = "libcurl")

options(download.file.method.GEOquery = "libcurl")


# ------------------------------------------------------------
# 4. DEFINE DATASET LOCATION
# ------------------------------------------------------------

geo_file <- "data/raw/GSE68801_series_matrix.txt.gz"

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE68nnn/GSE68801/matrix/",
  "GSE68801_series_matrix.txt.gz"
)


# ------------------------------------------------------------
# 5. DOWNLOAD GEO DATASET IF NOT ALREADY AVAILABLE
# ------------------------------------------------------------

if (!file.exists(geo_file)) {
  
  cat("\nDownloading GSE68801 dataset...\n")
  
  download.file(
    url = geo_url,
    destfile = geo_file,
    mode = "wb",
    method = "libcurl"
  )
  
} else {
  
  cat("\nPreviously downloaded GEO file found.\n")
  cat("Reusing existing file.\n")
  
}


# Check downloaded file

if (!file.exists(geo_file) ||
    file.info(geo_file)$size <= 0) {
  
  stop("GEO dataset download failed. Please check your internet connection.")
  
}

cat("\nGEO file successfully available.\n")

cat(
  "File size:",
  round(file.info(geo_file)$size / 1024^2, 2),
  "MB\n"
)


# ------------------------------------------------------------
# 6. IMPORT GEO DATASET
# ------------------------------------------------------------

cat("\nImporting GEO dataset...\n")

gse <- getGEO(
  filename = geo_file,
  GSEMatrix = TRUE,
  getGPL = FALSE
)


# ------------------------------------------------------------
# 7. EXTRACT EXPRESSIONSET OBJECT
# ------------------------------------------------------------

# The downloaded file is imported as an ExpressionSet.
# Do not use gse[[1]] in this case.

if (methods::is(gse, "ExpressionSet")) {
  
  eset <- gse
  
} else if (is.list(gse) && length(gse) > 0 &&
           methods::is(gse[[1]], "ExpressionSet")) {
  
  eset <- gse[[1]]
  
} else {
  
  stop("Unable to identify a valid ExpressionSet object.")
  
}


# ------------------------------------------------------------
# 8. EXTRACT GENE EXPRESSION MATRIX
# ------------------------------------------------------------

expression_data <- exprs(eset)

cat("\nExpression matrix extracted successfully.\n")


# ------------------------------------------------------------
# 9. EXTRACT SAMPLE METADATA
# ------------------------------------------------------------

sample_info <- pData(eset)

cat("Sample metadata extracted successfully.\n")


# ------------------------------------------------------------
# 10. CHECK DATASET DIMENSIONS
# ------------------------------------------------------------

cat("\n========== DATASET DIMENSIONS ==========\n")

cat("Number of probes:",
    nrow(expression_data), "\n")

cat("Number of samples:",
    ncol(expression_data), "\n")

cat("Number of metadata records:",
    nrow(sample_info), "\n")


# Verify expression and metadata sample counts

if (ncol(expression_data) != nrow(sample_info)) {
  
  stop("Expression matrix and sample metadata dimensions do not match.")
  
}


# ------------------------------------------------------------
# 11. DISPLAY SAMPLE INFORMATION
# ------------------------------------------------------------

cat("\n========== SAMPLE INFORMATION ==========\n")

cat("\nAvailable metadata columns:\n")

print(colnames(sample_info))


# Display sample titles

cat("\nSample titles:\n")

print(head(sample_info$title, 10))


# ------------------------------------------------------------
# 12. IDENTIFY SAMPLE CHARACTERISTIC COLUMNS
# ------------------------------------------------------------

characteristic_columns <- grep(
  "characteristics_ch1",
  colnames(sample_info),
  value = TRUE,
  ignore.case = TRUE
)


# Display sample characteristics

cat("\n========== SAMPLE CHARACTERISTICS ==========\n")

if (length(characteristic_columns) > 0) {
  
  print(
    head(
      sample_info[
        ,
        characteristic_columns,
        drop = FALSE
      ],
      10
    )
  )
  
} else {
  
  cat("No sample characteristic columns found.\n")
  
}


# ------------------------------------------------------------
# 13. SAVE EXPRESSION MATRIX
# ------------------------------------------------------------

saveRDS(
  expression_data,
  file = "data/raw/GSE68801_expression.rds"
)

cat("\nExpression matrix saved successfully.\n")


# ------------------------------------------------------------
# 14. SAVE SAMPLE METADATA
# ------------------------------------------------------------

write.csv(
  sample_info,
  file = "data/raw/GSE68801_sample_metadata.csv",
  row.names = TRUE
)

cat("Sample metadata saved successfully.\n")


# ------------------------------------------------------------
# 15. CREATE DATASET SUMMARY
# ------------------------------------------------------------

summary_text <- c(
  
  "ALOPECIA AREATA BIOINFORMATICS PROJECT",
  
  "---------------------------------------",
  
  "Dataset: GSE68801",
  
  "Source: NCBI Gene Expression Omnibus",
  
  paste(
    "Number of probes:",
    nrow(expression_data)
  ),
  
  paste(
    "Number of samples:",
    ncol(expression_data)
  ),
  
  paste(
    "Platform:",
    annotation(eset)
  ),
  
  paste(
    "Number of metadata records:",
    nrow(sample_info)
  )
  
)


writeLines(
  summary_text,
  "data/raw/dataset_summary.txt"
)


# ------------------------------------------------------------
# 16. VERIFY SAVED FILES
# ------------------------------------------------------------

cat("\n========== SAVED FILES ==========\n")

saved_files <- c(
  "data/raw/GSE68801_series_matrix.txt.gz",
  "data/raw/GSE68801_expression.rds",
  "data/raw/GSE68801_sample_metadata.csv",
  "data/raw/dataset_summary.txt"
)

for (file in saved_files) {
  
  cat(
    file,
    ":",
    file.exists(file),
    "\n"
  )
  
}


# ------------------------------------------------------------
# 17. FINAL COMPLETION MESSAGE
# ------------------------------------------------------------

cat("\n========================================\n")

cat("STAGE 1 DATA IMPORT COMPLETED!\n")

cat("Dataset: GSE68801\n")

cat("Probes:", nrow(expression_data), "\n")

cat("Samples:", ncol(expression_data), "\n")

cat("========================================\n")

#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages(library(tximport))

# --- Get Snakemake variables ---
quant_files <- snakemake@input$quant_files
tx2gene_path <- snakemake@input$tx2gene
counts_out <- snakemake@output$counts
tpm_out <- snakemake@output$tpm
sample_names <- snakemake@params$sample_names

# Set the names of the files vector to the sample names
names(quant_files) <- sample_names

# --- Read tx2gene map ---
cat("Loading tx2gene map...\n")
tx2gene <- read.delim(
  tx2gene_path,
  header = FALSE,
  col.names = c("TXNAME", "GENEID"),
  stringsAsFactors = FALSE,
  sep = "\t"
)
cat("Loaded transcript-to-gene map. First 5 rows:\n")
print(head(tx2gene))

# --- Run tximport ---
# Read in transcript-level counts and summarize to gene level
cat("Running tximport...\n")
txi <- tximport(quant_files,
                type = "kallisto",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

# --- Save output files ---
cat("tximport complete. Writing output matrices.\n")

# --- Write the raw integer counts ---
counts_df <- as.data.frame(round(txi$counts))
counts_matrix <- data.frame(gene_id = rownames(counts_df), counts_df)
write.table(
  counts_matrix,
  file = counts_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# --- Write the TPM matrix --
tpm_df <- as.data.frame(txi$abundance)
tpm_matrix <- data.frame(gene_id = rownames(tpm_df), tpm_df)
# Replaced write_tsv() with base R's write.table()
write.table(
  tpm_matrix,
  file = tpm_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("All done.\n")


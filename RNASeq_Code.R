# Install packages if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

library(biomaRt)
library(readxl)
library(writexl)

# File path
file_path <- "D:/RNASeq_project/GSE50499_raw_counts.xlsx"

# Read Excel file
counts_df <- read_excel(file_path)

# Connect to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map NCBI Gene IDs to Ensembl IDs
mapping <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  filters = "entrezgene_id",
  values = counts_df$GeneID,
  mart = ensembl
)

# Merge mapping into your counts table
counts_df <- merge(mapping, counts_df, by.x = "entrezgene_id", by.y = "GeneID", all.y = TRUE)

# Replace GeneID with Ensembl ID
counts_df <- counts_df[, c("ensembl_gene_id", colnames(counts_df)[3:ncol(counts_df)])]
colnames(counts_df)[1] <- "GeneID"

# Save updated file
write_xlsx(counts_df, "D:/RNASeq_project/GSE50499_raw_counts_with_ensembl.xlsx")

cat("✅ New file saved: GSE50499_raw_counts_with_ensembl.xlsx\n")






library(readxl)
library(writexl)

# Path to your Excel file
file_path <- "D:/RNASeq_project/counts_with_ensembl.xlsx"

# Read Excel file
df <- read_excel(file_path)

# Keep Gene IDs
gene_ids <- df$Geneid

# Clean column names: remove full paths and file extensions
colnames(df) <- sub(".*/", "", colnames(df))        # keep only file name
colnames(df) <- sub("_sorted.bam", "", colnames(df)) # remove suffix

# Extract count data only
count_data <- df[, -1]

# List of SRR IDs (two are technical replicates for each sample)
srr_ids <- c(
  "SRR960455", "SRR960456",
  "SRR960457", "SRR960458",
  "SRR960459", "SRR960460",
  "SRR960461", "SRR960462",
  "SRR960463", "SRR960464",
  "SRR960465", "SRR960466",
  "SRR960467", "SRR960468",
  "SRR960469", "SRR960470"
)

# Create sample names (Sample1, Sample2, etc.)
sample_names <- paste0("Sample", seq(1, length(srr_ids) / 2))

# Sum technical replicates
summed_counts <- data.frame(Geneid = gene_ids)

for (i in seq(1, length(srr_ids), by = 2)) {
  rep1 <- srr_ids[i]
  rep2 <- srr_ids[i + 1]
  sample_name <- paste0("Sample", (i + 1) / 2)
  
  summed_counts[[sample_name]] <- rowSums(count_data[, c(rep1, rep2)])
}

# Save new file
write_xlsx(summed_counts, "D:/RNASeq_project/counts_summed.xlsx")

cat("✅ Summed counts saved to: D:/RNASeq_project/counts_summed.xlsx\n")





# =========================
# Install & Load Packages
# =========================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}




library(DESeq2)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)

# =========================
# 1. Read Counts from Excel
# =========================
counts_df <- read_excel("D:/RNASeq_project/counts_summed.xlsx")

counts <- as.data.frame(counts_df)
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
counts <- round(counts)

# =========================
# 2. Create Sample Information
# =========================
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("MOV10_KD", "MOV10_KD",
                "MOV10_OE", "MOV10_OE", "MOV10_OE",
                "Control", "Control", "Control")
)
coldata$condition <- factor(coldata$condition,
                            levels = c("Control", "MOV10_KD", "MOV10_OE"))
stopifnot(all(colnames(counts) == rownames(coldata)))

# =========================
# 3. Create DESeq2 Dataset
# =========================
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# =========================
# 4. Get Shrunken Results
# =========================
res_KD <- lfcShrink(dds, coef = "condition_MOV10_KD_vs_Control", type = "apeglm") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

res_OE <- lfcShrink(dds, coef = "condition_MOV10_OE_vs_Control", type = "apeglm") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

# =========================
# 5. Clean & Classify
# =========================
clean_results <- function(res, padj_cut = 0.05, lfc_cut = 1) {
  res <- res[!is.na(res$log2FoldChange), ]
  res$padj[is.na(res$padj)] <- 1
  
  min_nonzero <- min(res$padj[res$padj > 0], na.rm = TRUE)
  if (is.finite(min_nonzero)) {
    res$padj[res$padj == 0] <- min_nonzero / 10
  }
  
  res$minusLog10Padj <- -log10(res$padj)
  res$minusLog10PadjCapped <- pmin(res$minusLog10Padj, 50)
  
  res$Significance <- ifelse(
    res$padj < padj_cut & abs(res$log2FoldChange) >= lfc_cut,
    "Significant", "Not significant"
  )
  res$Significance <- factor(res$Significance,
                             levels = c("Not significant", "Significant"))
  
  return(res)
}

res_KD <- clean_results(res_KD)
res_OE <- clean_results(res_OE)

library(EnhancedVolcano)

EnhancedVolcano(
  res_KD,
  lab = res_KD$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'MOV10 Knockdown vs Control',
  subtitle = 'Differential Expression',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 4,
  colAlpha = 0.8,
  col = c('grey70', 'skyblue', 'salmon', 'red3'),
  selectLab = head(res_KD$gene[order(res_KD$padj)], 10), # Top 10 most significant
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  boxedLabels = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  border = 'full',
  axisLabSize = 14,
  titleLabSize = 18
)
EnhancedVolcano(
  res_OE,
  lab = res_OE$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'MOV10 Overexpression vs Control',
  subtitle = 'Differential Expression',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 4,
  colAlpha = 0.8,
  col = c('purple', 'blue', 'orange', 'green'), # Four distinct colors
  selectLab = head(res_OE$gene[order(res_OE$padj)], 10),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  boxedLabels = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  border = 'full',
  axisLabSize = 14,
  titleLabSize = 18
)





library(DESeq2)
library(pheatmap)
library(dplyr)

# -----------------------------
# 1. Select significant genes
# -----------------------------
# Criteria: padj < 0.05 & |log2FC| >= 1
sig_KD <- res_KD %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  arrange(padj)

sig_OE <- res_OE %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  arrange(padj)

# Option 1: union of KD & OE genes
sig_genes <- union(sig_KD$gene, sig_OE$gene)

# -----------------------------
# 2. Variance Stabilizing Transformation
# -----------------------------
vsd <- vst(dds, blind = FALSE)

# Extract transformed counts for significant genes
mat <- assay(vsd)[sig_genes, ]

# -----------------------------
# 3. Scale genes for heatmap
# -----------------------------
mat_scaled <- t(scale(t(mat)))  # Z-score per gene

library(pheatmap)

annotation_col$Condition <- gsub("MOV10_KD", "KD", annotation_col$Condition)
annotation_col$Condition <- gsub("MOV10_OE", "OE", annotation_col$Condition)

# Re-run pheatmap
p <- pheatmap(mat_scaled,
              annotation_col = annotation_col,
              show_rownames = FALSE,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              fontsize = 10,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              main = "Expression Heatmap of Significant Genes")

grid.newpage()
grid.draw(p$gtable)



# =========================
# Install if needed
# =========================
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}

library(VennDiagram)
library(UpSetR)

# =========================
# Prepare DEG sets
# =========================
deg_list <- list(
  MOV10_KD = sig_KD$gene,
  MOV10_OE = sig_OE$gene
)


install.packages("ggVennDiagram")
library(ggVennDiagram)

ggVennDiagram(deg_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "skyblue") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(20, 50, 20, 100)
  ) +
  labs(title = "Overlap of DEGs: MOV10 KD vs OE") +
  coord_fixed(ratio = 1.0, clip = "off")   # makes circles narrower



# =========================
# 2. UpSetR Plot
# =========================
# Convert to binary membership table
upset_data <- fromList(deg_list)
UpSetR::upset(
  upset_data,
  nsets = length(deg_list),
  nintersects = NA,
  order.by = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "gray40",
  mainbar.y.label = "Number of Genes",
  sets.x.label = "Genes per Set"
)



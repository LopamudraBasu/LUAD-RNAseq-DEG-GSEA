# ==========================================================
# LUAD vs Normal Lung Bulk RNA-seq Analysis
# Differential Expression, GO, KEGG and GSEA
# Author: Lopamudra Basu
# ==========================================================

# ----------------------------------------------------------
# 1. Load Required Packages
# ----------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
# ----------------------------------------------------------
# 2. Define Project Directories
# ----------------------------------------------------------

project_dir <- getwd()

count_dir <- file.path(
  project_dir,
  "data",
  "counts"
)

results_dir <- file.path(
  project_dir,
  "results"
)

count_matrix_dir <- file.path(
  results_dir,
  "count_matrix"
)

dir.create(
  count_matrix_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

figures_dir <- file.path(
  results_dir,
  "figures"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ----------------------------------------------------------
# 3. Define Samples
# ----------------------------------------------------------

samples <- c(
  "SRR10440857",
  "SRR10440859",
  "SRR10440860",
  "SRR10440861",
  "SRR11357867",
  "SRR11357868",
  "SRR11357869",
  "SRR11357870"
)

# ----------------------------------------------------------
# 4. Import FeatureCounts Results
# ----------------------------------------------------------

count_list <- lapply(samples, function(s) {
  
  file <- file.path(
    count_dir,
    paste0(s, "_gene_counts.txt")
  )
  
  x <- read.delim(
    file,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  x[, c("Geneid", names(x)[ncol(x)])]
})

names(count_list) <- samples

for (i in seq_along(count_list)) {
  names(count_list[[i]])[2] <- samples[i]
}

# ----------------------------------------------------------
# 5. Combine Count Tables
# ----------------------------------------------------------

count_matrix <- Reduce(
  function(x, y)
    merge(
      x,
      y,
      by = "Geneid",
      all = FALSE
    ),
  count_list
)

dim(count_matrix)
head(count_matrix)
colnames(count_matrix)

# Save raw count matrix

write.table(
  count_matrix,
  file.path(
    count_matrix_dir,
    "raw_counts_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 6. Create Sample Metadata
# ----------------------------------------------------------

coldata <- data.frame(
  condition = c(
    "LUAD",
    "LUAD",
    "LUAD",
    "LUAD",
    "Normal",
    "Normal",
    "Normal",
    "Normal"
  ),
  row.names = samples
)

# Set Normal as reference condition

coldata$condition <- factor(
  coldata$condition,
  levels = c("Normal", "LUAD")
)

# Save metadata

write.csv(
  coldata,
  file.path(
    count_matrix_dir,
    "sample_metadata.csv"
  )
)

# ----------------------------------------------------------
# 7. Verify Count Matrix and Metadata
# ----------------------------------------------------------

stopifnot(
  all(
    colnames(count_matrix)[-1] ==
      rownames(coldata)
  )
)

# ----------------------------------------------------------
# 8. Prepare Count Matrix for DESeq2
# ----------------------------------------------------------

counts <- count_matrix[, -1]

rownames(counts) <- count_matrix$Geneid

# Basic data checks

stopifnot(!any(is.na(counts)))
stopifnot(!any(counts < 0))
stopifnot(all(counts == floor(counts)))

# ----------------------------------------------------------
# 9. Filter Lowly Expressed Genes
# ----------------------------------------------------------

keep <- rowSums(counts >= 10) >= 4

table(keep)

counts_filtered <- counts[keep, ]

dim(counts_filtered)

# ----------------------------------------------------------
# 10. Construct DESeq2 Dataset
# ----------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ condition
)

# ----------------------------------------------------------
# 11. Differential Expression Analysis
# ----------------------------------------------------------

dds <- DESeq(dds)

resultsNames(dds)

# ----------------------------------------------------------
# 12. Normalized Counts
# ----------------------------------------------------------

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

# ----------------------------------------------------------
# 13. Variance Stabilizing Transformation
# ----------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

# ----------------------------------------------------------
# 14. PCA Data
# ----------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(
    pca_data,
    "percentVar"
  )
)

percent_var

# ==========================================================
# PCA Visualization
# ==========================================================

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = name),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0.08)
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.12)
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of LUAD and Normal Lung Samples",
    x = paste0(
      "PC1: ",
      percent_var[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percent_var[2],
      "% variance"
    )
  )

ggsave(
  file.path(
    results_dir,
    "figures",
    "PCA_LUAD_vs_Normal.png"
  ),
  pca_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ==========================================================
# Differential Expression Analysis
# ==========================================================

res <- results(
  dds,
  contrast = c(
    "condition",
    "LUAD",
    "Normal"
  ),
  alpha = 0.05
)

res_ordered <- res[
  order(res$padj),
]

res_df <- as.data.frame(
  res_ordered
)

res_df$Geneid <- rownames(res_df)

# Remove Ensembl version suffix
res_df$ENSEMBL <- sub(
  "\\..*$",
  "",
  res_df$Geneid
)

# ==========================================================
# DEG Classification
# ==========================================================

deg_sig <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05,
]

deg_fc2 <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
]

deg_up <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange >= 1,
]

deg_down <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange <= -1,
]

# ==========================================================
# Save DEG Tables
# ==========================================================

write.csv(
  res_df,
  file.path(
    results_dir,
    "DEG_LUAD_vs_Normal_all_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_sig,
  file.path(
    results_dir,
    "DEG_LUAD_vs_Normal_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_fc2,
  file.path(
    results_dir,
    "DEG_LUAD_vs_Normal_FC2.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_up,
  file.path(
    results_dir,
    "DEG_LUAD_vs_Normal_upregulated.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_down,
  file.path(
    results_dir,
    "DEG_LUAD_vs_Normal_downregulated.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# Volcano Plot
# ==========================================================

volcano_df <- as.data.frame(res)

volcano_df$Geneid <- rownames(
  volcano_df
)

volcano_df <- volcano_df[
  !is.na(volcano_df$padj),
]

volcano_df$category <- "Not significant"

volcano_df$category[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange >= 1
] <- "Upregulated"

volcano_df$category[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange <= -1
] <- "Downregulated"

volcano_df$neg_log10_padj <- -log10(
  volcano_df$padj
)

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = category
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    title = "Differential Gene Expression: LUAD vs Normal Lung",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Category"
  )

ggsave(
  file.path(
    results_dir,
    "figures",
    "Volcano_LUAD_vs_Normal.png"
  ),
  volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)


# ============================================================
# Gene Annotation and Top 50 DEG Heatmap
# ============================================================

# ------------------------------------------------------------
# 1. Clean Ensembl IDs
# ------------------------------------------------------------

res_df$ENSEMBL <- sub("\\..*$", "", res_df$Geneid)

# ------------------------------------------------------------
# 2. Annotate Ensembl IDs with gene symbols
# ------------------------------------------------------------

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res_df$ENSEMBL),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

# Keep one symbol per Ensembl ID
gene_map <- gene_map[
  !duplicated(gene_map$ENSEMBL),
]

res_annotated <- merge(
  res_df,
  gene_map,
  by = "ENSEMBL",
  all.x = TRUE
)

write.csv(
  res_annotated,
  "results/DEG_LUAD_vs_Normal_annotated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Select top 50 genes based on adjusted p-value
# ------------------------------------------------------------

res_valid <- res_df[
  !is.na(res_df$padj),
]

top50 <- res_valid[
  order(res_valid$padj),
][1:50, ]

# Add gene symbols
top50$SYMBOL <- res_annotated$SYMBOL[
  match(
    top50$ENSEMBL,
    res_annotated$ENSEMBL
  )
]

# Check annotation
sum(is.na(top50$SYMBOL))

# ------------------------------------------------------------
# 4. Extract VST expression matrix
# ------------------------------------------------------------

vsd_mat <- assay(vsd)

heatmap_mat <- vsd_mat[top50$Geneid, ]

# Use gene symbols as row names
rownames(heatmap_mat) <- top50$SYMBOL

# ------------------------------------------------------------
# 5. Row-wise Z-score scaling
# ------------------------------------------------------------

heatmap_mat_scaled <- t(
  scale(t(heatmap_mat))
)

# ------------------------------------------------------------
# 6. Sample annotation
# ------------------------------------------------------------

annotation_col <- data.frame(
  Condition = coldata$condition
)

rownames(annotation_col) <- rownames(coldata)

# ------------------------------------------------------------
# 7. Generate Top 50 DEG heatmap
# ------------------------------------------------------------

png(
  filename = "results/figures/Top50_DEG_Heatmap_LUAD_vs_Normal.png",
  width = 2400,
  height = 2800,
  res = 300
)

pheatmap(
  heatmap_mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 11,
  main = "Top 50 Differentially Expressed Genes\nLUAD vs Normal Lung",
  border_color = NA
)

dev.off()

# ============================================================
# Functional Enrichment Analysis
# GO Biological Process + KEGG
# ============================================================

# ------------------------------------------------------------
# 1. Define upregulated and downregulated genes
# ------------------------------------------------------------

up_genes <- res_annotated$ENSEMBL[
  !is.na(res_annotated$padj) &
    res_annotated$padj < 0.05 &
    res_annotated$log2FoldChange >= 1
]

down_genes <- res_annotated$ENSEMBL[
  !is.na(res_annotated$padj) &
    res_annotated$padj < 0.05 &
    res_annotated$log2FoldChange <= -1
]

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")


# ------------------------------------------------------------
# 2. Convert Ensembl IDs to Entrez IDs
# ------------------------------------------------------------

up_entrez <- bitr(
  unique(up_genes),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_entrez <- bitr(
  unique(down_genes),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat("Mapped upregulated genes:", nrow(up_entrez), "\n")
cat("Mapped downregulated genes:", nrow(down_entrez), "\n")


# ------------------------------------------------------------
# 3. Define background gene universe
# ------------------------------------------------------------

background_genes <- sub(
  "\\..*$",
  "",
  rownames(dds)
)

background_entrez <- bitr(
  unique(background_genes),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_ids <- unique(up_entrez$ENTREZID)
down_ids <- unique(down_entrez$ENTREZID)
background_ids <- unique(background_entrez$ENTREZID)

cat("Unique upregulated Entrez IDs:", length(up_ids), "\n")
cat("Unique downregulated Entrez IDs:", length(down_ids), "\n")
cat("Background Entrez IDs:", length(background_ids), "\n")


# ------------------------------------------------------------
# 4. GO Biological Process enrichment
# ------------------------------------------------------------

ego_up_BP <- enrichGO(
  gene = up_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_down_BP <- enrichGO(
  gene = down_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# ------------------------------------------------------------
# 5. Save GO enrichment results
# ------------------------------------------------------------

write.csv(
  as.data.frame(ego_up_BP),
  "results/GO_BP_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down_BP),
  "results/GO_BP_Downregulated.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. GO enrichment plots
# ------------------------------------------------------------

png(
  "results/figures/GO_BP_Upregulated_Top20.png",
  width = 2200,
  height = 1800,
  res = 250
)

dotplot(
  ego_up_BP,
  showCategory = 20,
  title = "GO Biological Process Enrichment - Upregulated Genes"
)

dev.off()


png(
  "results/figures/GO_BP_Downregulated_Top15.png",
  width = 2400,
  height = 1900,
  res = 250
)

dotplot(
  ego_down_BP,
  showCategory = 15,
  title = "GO Biological Process Enrichment - Downregulated Genes"
)

dev.off()


# ------------------------------------------------------------
# 7. KEGG enrichment
# ------------------------------------------------------------

ekegg_up <- enrichKEGG(
  gene = up_ids,
  universe = background_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

ekegg_down <- enrichKEGG(
  gene = down_ids,
  universe = background_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)


# ------------------------------------------------------------
# 8. Save KEGG enrichment results
# ------------------------------------------------------------

write.csv(
  as.data.frame(ekegg_up),
  "results/KEGG_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg_down),
  "results/KEGG_Downregulated.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. KEGG dotplots
# ------------------------------------------------------------

p_kegg_up <- dotplot(
  ekegg_up,
  showCategory = 7,
  title = "KEGG Pathway Enrichment - Upregulated Genes"
) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  "results/figures/KEGG_Upregulated_Dotplot.png",
  p_kegg_up,
  width = 10,
  height = 7,
  dpi = 300
)


p_kegg_down <- dotplot(
  ekegg_down,
  showCategory = 15,
  title = "KEGG Pathway Enrichment - Downregulated Genes"
) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  "results/figures/KEGG_Downregulated_Dotplot.png",
  p_kegg_down,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 10. Make KEGG results readable with gene symbols
# ------------------------------------------------------------

ekegg_up_readable <- setReadable(
  ekegg_up,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

ekegg_down_readable <- setReadable(
  ekegg_down,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# ============================================================
# KEGG Pathway-Gene Network Analysis
# ============================================================

# ------------------------------------------------------------
# 1. Prepare upregulated fold-change vector
# ------------------------------------------------------------

up_df <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange >= 1,
]

up_fc_map <- bitr(
  up_df$ENSEMBL,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_fc_map <- up_fc_map[
  !duplicated(up_fc_map$ENSEMBL),
]

up_fc_final <- up_df$log2FoldChange[
  match(
    up_fc_map$ENSEMBL,
    up_df$ENSEMBL
  )
]

names(up_fc_final) <- up_fc_map$ENTREZID
# ------------------------------------------------------------
# 2. Upregulated KEGG cnetplot
# ------------------------------------------------------------

p_cnet_up <- cnetplot(
  ekegg_up_readable,
  showCategory = 5,
  foldChange = up_fc_final,
  layout = "fr",
  node_label = "all"
) +
  ggtitle(
    "KEGG Pathway-Gene Network - Upregulated Genes"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

ggsave(
  "results/figures/KEGG_cnetplot_upregulated.png",
  p_cnet_up,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------
# 3. Prepare downregulated fold-change vector
# ------------------------------------------------------------

down_df <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange <= -1,
]

down_fc_map <- bitr(
  down_df$ENSEMBL,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_fc_map <- down_fc_map[
  !duplicated(down_fc_map$ENSEMBL),
]

down_fc_final <- down_df$log2FoldChange[
  match(
    down_fc_map$ENSEMBL,
    down_df$ENSEMBL
  )
]

names(down_fc_final) <- down_fc_map$ENTREZID


# ------------------------------------------------------------
# 4. Downregulated KEGG cnetplot
# ------------------------------------------------------------

p_cnet_down <- cnetplot(
  ekegg_down_readable,
  showCategory = 5,
  foldChange = down_fc_final,
  layout = "fr",
  node_label = "all"
) +
  ggtitle(
    "KEGG Pathway-Gene Network - Downregulated Genes"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

ggsave(
  "results/figures/KEGG_cnetplot_downregulated.png",
  p_cnet_down,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)


# ============================================================
# GSEA Preparation
# ============================================================

# ------------------------------------------------------------
# 5. Create ranked gene list using DESeq2 Wald statistic
# ------------------------------------------------------------

gsea_df <- res_df[
  !is.na(res_df$stat) &
    !is.na(res_df$ENSEMBL),
]

gsea_df <- gsea_df[
  !duplicated(gsea_df$ENSEMBL),
]

gene_list <- gsea_df$stat

names(gene_list) <- gsea_df$ENSEMBL

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


# ------------------------------------------------------------
# 6. Convert Ensembl IDs to Entrez IDs
# ------------------------------------------------------------

gsea_map <- bitr(
  names(gene_list),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gsea_map <- gsea_map[
  !duplicated(gsea_map$ENSEMBL),
]


# ============================================================
# GSEA - KEGG Pathway Enrichment
# ============================================================

# ------------------------------------------------------------
# 1. Prepare ranked Ensembl gene list
# ------------------------------------------------------------

gsea_rank_df <- data.frame(
  ENSEMBL = names(gene_list),
  stat = as.numeric(gene_list),
  stringsAsFactors = FALSE
)

# Add Entrez IDs
gsea_rank_df <- merge(
  gsea_rank_df,
  gsea_map,
  by = "ENSEMBL"
)

# Keep the strongest statistic for duplicated Entrez IDs
gsea_rank_df <- gsea_rank_df[
  order(
    gsea_rank_df$ENTREZID,
    -abs(gsea_rank_df$stat)
  ),
]

gsea_rank_df <- gsea_rank_df[
  !duplicated(gsea_rank_df$ENTREZID),
]

# Create final ranked Entrez vector
gene_list_entrez <- gsea_rank_df$stat

names(gene_list_entrez) <- gsea_rank_df$ENTREZID

gene_list_entrez <- sort(
  gene_list_entrez,
  decreasing = TRUE
)

# Validate ranking
stopifnot(
  !anyDuplicated(names(gene_list_entrez))
)

stopifnot(
  !any(is.na(names(gene_list_entrez)))
)

cat(
  "Genes in GSEA ranking:",
  length(gene_list_entrez),
  "\n"
)


# ------------------------------------------------------------
# 2. Run KEGG GSEA
# ------------------------------------------------------------

gsea_kegg <- gseKEGG(
  geneList = gene_list_entrez,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)


# ------------------------------------------------------------
# 3. Convert GSEA results to data frame
# ------------------------------------------------------------

gsea_res <- as.data.frame(gsea_kegg)

gsea_res <- gsea_res[
  order(gsea_res$NES, decreasing = TRUE),
]

write.csv(
  gsea_res,
  "results/GSEA_KEGG_results.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. GSEA dotplot
# ------------------------------------------------------------

p_gsea_dot <- dotplot(
  gsea_kegg,
  showCategory = 20,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  ggtitle(
    "KEGG GSEA - LUAD vs Normal Lung"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

ggsave(
  "results/figures/KEGG_GSEA_dotplot.png",
  p_gsea_dot,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------
# 5. Top positively enriched pathways
# ------------------------------------------------------------

gsea_positive <- gsea_res[
  gsea_res$NES > 0,
]

gsea_positive <- gsea_positive[
  order(gsea_positive$NES, decreasing = TRUE),
]

gsea_positive <- head(
  gsea_positive,
  10
)

gsea_positive$Description <- factor(
  gsea_positive$Description,
  levels = rev(gsea_positive$Description)
)

p_gsea_positive <- ggplot(
  gsea_positive,
  aes(
    x = NES,
    y = Description,
    size = -log10(p.adjust),
    color = NES
  )
) +
  geom_point() +
  labs(
    title = "GSEA: Top Positively Enriched KEGG Pathways",
    subtitle = "LUAD vs Normal Lung",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    panel.grid.major.y = element_blank()
  )

ggsave(
  "results/figures/GSEA_KEGG_Positive_NES_Top10.png",
  p_gsea_positive,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------
# 6. Top negatively enriched pathways
# ------------------------------------------------------------

gsea_negative <- gsea_res[
  gsea_res$NES < 0,
]

gsea_negative <- gsea_negative[
  order(gsea_negative$NES, decreasing = FALSE),
]

gsea_negative <- head(
  gsea_negative,
  10
)

gsea_negative$Description <- factor(
  gsea_negative$Description,
  levels = rev(gsea_negative$Description)
)

p_gsea_negative <- ggplot(
  gsea_negative,
  aes(
    x = NES,
    y = Description,
    size = -log10(p.adjust),
    color = NES
  )
) +
  geom_point() +
  labs(
    title = "GSEA: Top Negatively Enriched KEGG Pathways",
    subtitle = "LUAD vs Normal Lung",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    panel.grid.major.y = element_blank()
  )

ggsave(
  "results/figures/GSEA_KEGG_Negative_NES_Top10.png",
  p_gsea_negative,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ============================================================
# GSEA Individual Pathway Analysis
# ============================================================

# ------------------------------------------------------------
# 1. Individual GSEA enrichment plots
# ------------------------------------------------------------

p_cellcycle <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04110",
  title = "GSEA - Cell Cycle",
  pvalue_table = TRUE
)

ggsave(
  "results/figures/GSEA_Cell_Cycle.png",
  p_cellcycle,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


p_dna <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa03030",
  title = "GSEA - DNA Replication",
  pvalue_table = TRUE
)

ggsave(
  "results/figures/GSEA_DNA_Replication.png",
  p_dna,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


p_fanconi <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa03460",
  title = "GSEA - Fanconi Anemia Pathway",
  pvalue_table = TRUE
)

ggsave(
  "results/figures/GSEA_Fanconi_Anemia.png",
  p_fanconi,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


p_cgmp <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04022",
  title = "GSEA - cGMP-PKG Signaling Pathway",
  pvalue_table = TRUE
)

ggsave(
  "results/figures/GSEA_cGMP_PKG.png",
  p_cgmp,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


p_hormone <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04081",
  title = "GSEA - Hormone Signaling",
  pvalue_table = TRUE
)

ggsave(
  "results/figures/GSEA_Hormone_Signaling.png",
  p_hormone,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


# ============================================================
# Leading-Edge / Core-Enrichment Analysis
# ============================================================

key_pathways <- c(
  "hsa04110",
  "hsa03030",
  "hsa03460",
  "hsa04022",
  "hsa04081"
)

gsea_key <- gsea_kegg@result[
  match(key_pathways, gsea_kegg@result$ID),
  c(
    "ID",
    "Description",
    "NES",
    "pvalue",
    "p.adjust",
    "core_enrichment"
  )
]

write.csv(
  gsea_key,
  "results/GSEA_key_pathways_leading_edge.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 2. Convert core enrichment Entrez IDs to gene symbols
# ------------------------------------------------------------

core_gene_tables <- lapply(
  seq_len(nrow(gsea_key)),
  function(i) {
    
    genes <- unlist(
      strsplit(
        gsea_key$core_enrichment[i],
        "/"
      )
    )
    
    gene_map <- bitr(
      genes,
      fromType = "ENTREZID",
      toType = "SYMBOL",
      OrgDb = org.Hs.eg.db
    )
    
    gene_map <- gene_map[
      !duplicated(gene_map$ENTREZID),
    ]
    
    gene_map$Pathway <- gsea_key$Description[i]
    
    gene_map
  }
)

core_genes_all <- do.call(
  rbind,
  core_gene_tables
)


# ------------------------------------------------------------
# 3. Save leading-edge gene membership
# ------------------------------------------------------------

core_gene_summary <- aggregate(
  SYMBOL ~ Pathway,
  data = core_genes_all,
  FUN = function(x) paste(x, collapse = "/")
)

write.csv(
  core_gene_summary,
  "results/GSEA_leading_edge_genes.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Identify genes shared between selected pathways
# ------------------------------------------------------------

core_gene_list <- split(
  core_genes_all$SYMBOL,
  core_genes_all$Pathway
)

gene_frequency <- sort(
  table(core_genes_all$SYMBOL),
  decreasing = TRUE
)

shared_core_genes <- gene_frequency[
  gene_frequency >= 2
]

shared_core_df <- data.frame(
  SYMBOL = names(shared_core_genes),
  Pathway_Count = as.integer(shared_core_genes)
)

write.csv(
  shared_core_df,
  "results/GSEA_shared_leading_edge_genes.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Pathway membership of shared genes
# ------------------------------------------------------------

shared_gene_pathways <- core_genes_all[
  core_genes_all$SYMBOL %in%
    names(shared_core_genes),
]

shared_gene_pathways <- shared_gene_pathways[
  order(
    shared_gene_pathways$SYMBOL,
    shared_gene_pathways$Pathway
  ),
]

shared_pathway_summary <- aggregate(
  Pathway ~ SYMBOL,
  data = shared_gene_pathways,
  FUN = paste,
  collapse = " | "
)

write.csv(
  shared_pathway_summary,
  "results/GSEA_shared_gene_pathway_membership.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Add DESeq2 statistics to shared leading-edge genes
# ------------------------------------------------------------

shared_le_genes <- merge(
  shared_core_df,
  shared_pathway_summary,
  by = "SYMBOL",
  all.x = TRUE
)

symbol_entrez <- unique(
  core_genes_all[, c(
    "SYMBOL",
    "ENTREZID"
  )]
)

shared_le_genes <- merge(
  shared_le_genes,
  symbol_entrez,
  by = "SYMBOL",
  all.x = TRUE
)

gene_stats <- unique(
  res_annotated[, c(
    "SYMBOL",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )]
)

shared_le_genes <- merge(
  shared_le_genes,
  gene_stats,
  by = "SYMBOL",
  all.x = TRUE
)

shared_le_genes <- shared_le_genes[
  order(
    abs(shared_le_genes$log2FoldChange),
    decreasing = TRUE
  ),
]

write.csv(
  shared_le_genes,
  "results/GSEA_shared_leading_edge_ranked_genes.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Candidate leading-edge genes
# ------------------------------------------------------------

candidate_genes <- shared_le_genes[
  order(
    shared_le_genes$padj,
    -abs(shared_le_genes$log2FoldChange)
  ),
]

write.csv(
  candidate_genes,
  "results/GSEA_leading_edge_candidate_genes.csv",
  row.names = FALSE
)

cellcycle_dna_overlap <- intersect(
  core_gene_list[["Cell cycle"]],
  core_gene_list[["DNA replication"]]
)

cellcycle_fanconi_overlap <- intersect(
  core_gene_list[["Cell cycle"]],
  core_gene_list[["Fanconi anemia pathway"]]
)

dna_fanconi_overlap <- intersect(
  core_gene_list[["DNA replication"]],
  core_gene_list[["Fanconi anemia pathway"]]
)

overlap_summary <- data.frame(
  Comparison = c(
    "Cell cycle vs DNA replication",
    "Cell cycle vs Fanconi anemia",
    "DNA replication vs Fanconi anemia"
  ),
  Shared_Genes = c(
    paste(cellcycle_dna_overlap, collapse = ", "),
    paste(cellcycle_fanconi_overlap, collapse = ", "),
    paste(dna_fanconi_overlap, collapse = ", ")
  )
)

write.csv(
  overlap_summary,
  "results/GSEA_pathway_overlap_summary.csv",
  row.names = FALSE
)
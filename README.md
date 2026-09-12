# LUAD Bulk RNA-seq Differential Expression and Pathway Analysis

## Overview

This project presents an end-to-end bulk RNA-seq analysis of **Lung Adenocarcinoma (LUAD)** versus **Normal lung tissue** using publicly available RNA-seq data from the NCBI Gene Expression Omnibus (GEO).

The analysis combines command-line NGS processing with downstream statistical and functional analysis in R.

The workflow includes:

- RNA-seq data processing
- Quality control and read trimming
- Reference genome alignment
- SAMtools processing
- Gene-level read quantification
- Count matrix construction
- Differential expression analysis using DESeq2
- Principal Component Analysis (PCA)
- Differentially expressed gene (DEG) identification
- Gene Ontology (GO) enrichment
- KEGG pathway enrichment
- KEGG pathway-gene network analysis
- KEGG Gene Set Enrichment Analysis (GSEA)
- Leading-edge/core-enrichment analysis
- Shared pathway gene analysis

---

## Biological Question

The primary objective was to identify genes and biological pathways that differ between **LUAD tumor tissue and normal lung tissue**.

The analysis was designed to answer:

1. Which genes are differentially expressed between LUAD and normal tissue?
2. Which biological processes and pathways are enriched among the differentially expressed genes?
3. Which pathways show coordinated changes across the complete ranked gene set?
4. Which leading-edge genes contribute to the enrichment of important pathways?

---

## Dataset

**GEO Accession:** GSE140343

The original study contains LUAD tumor samples and paired non-cancerous adjacent tissue samples.

For this analysis, **8 samples** were selected:

- 4 LUAD samples
- 4 Normal samples

### Selected Samples

| Condition | Sample IDs |
|---|---|
| LUAD | SRR10440857, SRR10440859, SRR10440860, SRR10440861 |
| Normal | SRR11357867, SRR11357868, SRR11357869, SRR11357870 |

The selected samples were analyzed as two experimental conditions: **LUAD** and **Normal**.

---

# Analysis Workflow

```text
RNA-seq data
     │
     ▼
FASTQ generation
     │
     ▼
Quality Control
     │
     ▼
Read trimming
     │
     ▼
Reference genome alignment
     │
     ▼
SAMtools processing
     │
     ▼
Gene-level quantification
     │
     ▼
Raw count matrix
     │
     ▼
Low-count filtering
     │
     ▼
DESeq2 analysis
     │
     ├── Normalization
     ├── Differential expression
     └── Statistical testing
     │
     ▼
Sample-level QC
     ├── PCA
     ├── Sample correlation
     └── Sample-to-sample distance
     │
     ▼
DEG identification
     │
     ├── Volcano plot
     └── Top DEG heatmap
     │
     ▼
Functional enrichment
     ├── GO Biological Process
     └── KEGG ORA
     │
     ▼
KEGG pathway-gene network analysis
     │
     ▼
KEGG GSEA
     │
     ▼
Leading-edge analysis
     │
     ▼
Shared pathway genes
1. RNA-seq Data Processing

The sequencing data were processed locally using command-line bioinformatics tools.

The processing workflow included:

SRA Toolkit for sequencing-data retrieval
FASTQ generation
Fastp for quality control and read trimming
HISAT2 for reference genome alignment
SAMtools for alignment processing
featureCounts for gene-level read quantification

The resulting gene-level count files were used as input for the downstream R analysis.

Large sequencing and reference files were intentionally not included in the GitHub repository.

2. Count Matrix Construction

Gene-level count files from the selected samples were imported into R and combined into a single count matrix.

The analysis used:

Gene identifiers as rows
Samples as columns
Raw integer read counts as input

The sample metadata were defined according to the two experimental conditions:

Normal
LUAD

Normal tissue was used as the reference condition.

The differential expression contrast was:

LUAD vs Normal
3. Count Filtering

Low-abundance genes were removed before differential expression analysis.

The filtering criterion was:

At least 10 counts in at least 4 samples

This reduces the influence of genes with extremely low read counts while retaining genes with sufficient evidence for downstream statistical analysis.

4. Differential Expression Analysis

Differential expression analysis was performed using DESeq2.

The experimental design was:

~ condition

DESeq2 was used for:

Library-size normalization
Estimation of dispersion
Model fitting
Wald statistical testing
Multiple-testing correction

The differential expression contrast was:

LUAD vs Normal
Significance Criteria

Genes were considered statistically significant when:

Adjusted p-value < 0.05

For fold-change-based DEG analysis, an additional threshold was applied:

|log2FoldChange| >= 1

Therefore:

Upregulated genes
log2FoldChange >= 1
padj < 0.05
Downregulated genes
log2FoldChange <= -1
padj < 0.05

Positive log2 fold-change values represent higher expression in LUAD relative to Normal, while negative values represent higher expression in Normal relative to LUAD.

5. Sample-Level Analysis and Visualization

Several visualization approaches were used to evaluate sample-level structure and expression patterns.

Principal Component Analysis

PCA was performed using variance-stabilizing transformed expression data.

PCA was used to assess:

Sample clustering
Separation between LUAD and Normal groups
Major sources of variation

The variance-stabilized data were used for visualization rather than as input for the DESeq2 statistical test.

Sample Correlation

Pairwise sample correlations were examined to evaluate the similarity of global expression profiles between samples.

Sample-to-Sample Distance

Sample-to-sample distances were calculated to visualize relationships among the eight samples.

DEG Heatmap

A heatmap of the top 50 genes ranked by adjusted p-value was generated using variance-stabilized expression values.

Row-wise scaling was applied for visualization.

Volcano Plot

A volcano plot was generated to visualize:

Log2 fold-change
Statistical significance
Upregulated genes
Downregulated genes
6. Gene Annotation

Ensembl gene identifiers were processed to remove version suffixes before annotation.

Gene identifiers were mapped to gene symbols using:

org.Hs.eg.db
AnnotationDbi

The annotated differential-expression results were saved for downstream interpretation.

7. Gene Ontology Enrichment

GO Biological Process enrichment analysis was performed separately for:

Upregulated genes
Downregulated genes

The enrichment analysis was performed using clusterProfiler.

The background universe consisted of genes retained after count filtering and successfully mapped to Entrez identifiers.

Multiple testing was controlled using the Benjamini-Hochberg (BH) method.

Top enriched biological processes were visualized using dot plots.

8. KEGG Over-Representation Analysis

KEGG pathway enrichment was performed separately for:

Upregulated genes
Downregulated genes

The analysis used:

clusterProfiler
org.Hs.eg.db

The background gene universe was derived from the genes retained in the filtered DESeq2 dataset.

KEGG pathway results were visualized using:

Dot plots
Pathway-gene network plots (cnetplots)

The pathway-gene networks provide a view of the relationships between enriched pathways and the genes contributing to those pathways.

9. KEGG Pathway-Gene Network Analysis

KEGG cnetplots were generated for the enriched upregulated and downregulated pathways.

The networks combine:

Enriched pathways
Contributing genes
Gene-level fold-change information

This provides a more detailed view of how individual genes contribute to multiple enriched pathways.

10. Gene Set Enrichment Analysis (GSEA)

KEGG GSEA was performed using the complete ranked gene list rather than restricting the analysis to statistically significant DEGs.

Genes were ranked using the DESeq2 Wald statistic.

This approach allows pathway-level patterns to be detected even when individual genes do not meet the predefined DEG threshold.

The ranked gene list was mapped to Entrez identifiers before GSEA.

The ranking was checked to ensure that the final Entrez-based gene list contained unique identifiers.

11. GSEA Pathway Analysis

Several biologically relevant KEGG pathways were investigated individually using GSEA enrichment plots.

Selected pathways included:

Cell cycle
DNA replication
Fanconi anemia pathway
cGMP-PKG signaling pathway
Hormone signaling

Both positively and negatively enriched pathways were examined based on their normalized enrichment scores (NES).

A positive NES indicates enrichment toward the LUAD end of the ranked gene list, while a negative NES indicates enrichment toward the Normal end.

12. Leading-Edge Analysis

Leading-edge/core-enrichment analysis was performed on selected GSEA pathways.

The analysis extracted the genes contributing most strongly to the enrichment signal of each selected pathway.

The workflow included:

Extraction of GSEA core-enrichment genes
Conversion of Entrez identifiers to gene symbols
Identification of genes shared between pathways
Counting pathway membership
Integration with differential-expression statistics
Ranking candidate shared genes

This provides an additional layer of interpretation beyond pathway-level enrichment.

Importantly, a gene can contribute to GSEA enrichment even if it does not independently satisfy the predefined DEG cutoff.

13. Key Findings

The analysis identified substantial transcriptional differences between LUAD and Normal lung tissue.

Several pathways were positively enriched toward the LUAD side of the ranked gene list, including pathways associated with:

Cell cycle
DNA replication
Fanconi anemia
Mismatch repair
Homologous recombination
Base excision repair
Ribosome biogenesis
Nucleotide metabolism

Several signaling-related pathways showed negative enrichment, including:

cGMP-PKG signaling
Hormone signaling
Calcium signaling
Adrenergic signaling
Vascular smooth muscle contraction
Renin secretion

The enrichment results indicate strong differences in proliferative and DNA-repair-associated biological programs between LUAD and Normal tissue, together with changes in signaling-related pathways.

14. Shared Leading-Edge Genes

The leading-edge analysis identified genes shared across multiple selected pathways.

Examples of genes shared between the Cell Cycle and DNA Replication pathways included:

MCM2
MCM3
MCM4
MCM5
MCM6
PCNA

The analysis also identified shared genes between the cGMP-PKG and Hormone signaling pathways.

These shared genes were further integrated with differential-expression statistics to prioritize candidate genes contributing to multiple pathway-level signals.

Repository Structure
LUAD-RNAseq-DEG-GSEA/
│
├── results/
│   ├── count_matrix/
│   ├── figures/
│   └── tables/
│
├── scripts/
│   └── LUAD_DESeq2_GO_KEGG_GSEA.R
│
├── README.md
└── LICENSE
Results Included

The repository contains the major outputs generated during the analysis, including:

Differential Expression
Complete DESeq2 results
Significant DEG table
FC2 DEG table
Upregulated DEG table
Downregulated DEG table
Annotated DEG table
Quality Control and Visualization
PCA plot
Sample correlation plot
Sample-to-sample distance plot
Top-50 DEG heatmap
Volcano plot
Functional Enrichment
GO Biological Process results
KEGG enrichment results
KEGG enrichment dot plots
KEGG pathway-gene network plots
GSEA
KEGG GSEA results
GSEA dot plot
Positive and negative NES plots
Individual pathway enrichment plots
Leading-Edge Analysis
Key pathway leading-edge results
Leading-edge gene tables
Shared leading-edge genes
Shared gene-pathway membership
Ranked shared leading-edge genes
Candidate leading-edge genes
Reproducibility

The complete downstream R analysis is provided in:

scripts/LUAD_DESeq2_GO_KEGG_GSEA.R

The script was tested from a fresh R session using the project root as the working directory.

The complete downstream workflow successfully executed from:

Count matrix construction
        ↓
DESeq2
        ↓
PCA
        ↓
DEG analysis
        ↓
GO enrichment
        ↓
KEGG enrichment
        ↓
KEGG cnetplot
        ↓
GSEA
        ↓
Individual GSEA pathway analysis
        ↓
Leading-edge analysis
        ↓
Shared leading-edge gene analysis

Large raw sequencing files, alignment files, genome references, and other computationally large files are not included in the repository.

Software and Tools
Command-Line Tools
SRA Toolkit
Fastp
HISAT2
SAMtools
featureCounts
R / Bioconductor
DESeq2
clusterProfiler
enrichplot
org.Hs.eg.db
AnnotationDbi
SummarizedExperiment
R Visualization Packages
ggplot2
ggrepel
pheatmap
Skills Demonstrated

This project demonstrates practical experience in:

Bulk RNA-seq analysis
NGS data processing
Linux command-line workflows
Quality control and read preprocessing
Genome alignment
Read quantification
Count-matrix construction
Experimental design
DESeq2 differential expression analysis
Multiple-testing correction
RNA-seq visualization
Gene annotation
GO enrichment
KEGG enrichment
Network-based pathway visualization
GSEA
Leading-edge analysis
Reproducible R scripting
Biological interpretation of transcriptomic data
Author

Lopamudra Basu

Bioinformatics Research Associate | Computational Genomics | RNA-seq | NGS

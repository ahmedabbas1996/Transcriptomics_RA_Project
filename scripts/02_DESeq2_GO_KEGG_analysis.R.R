# ============================================================
# Transcriptomics RA - Final script
# Startpunt: officiële count matrix (count_matrix_RA.txt)
# Auteur: Ahmed Abbas
#
# Dit script bevat alleen de stappen en figuren die in de GitHub-pagina
# worden gebruikt:
# 1. Count matrix inladen
# 2. Sample metadata maken
# 3. DESeq2 analyse
# 4. PCA plot
# 5. Volcano plot
# 6. GO analyse
# 7. KEGG analyse
# 8. Pathview voor TNF signaling pathway
# ============================================================


# ----------------------------
# 0. Clean environment
# ----------------------------

rm(list = ls())


# ----------------------------
# 1. Packages installeren/laden
# ----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

options(timeout = 10000)
options(BioC_mirror = "https://bioconductor.statistik.tu-dortmund.de")

bioc_packages <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot",
  "pathview"
)

cran_packages <- c(
  "ggplot2",
  "ggrepel"
)

for (p in bioc_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

for (p in cran_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(pathview)


# ----------------------------
# 2. Mappen aanmaken
# ----------------------------

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ----------------------------
# 3. Count matrix inladen
# ----------------------------

counts <- read.table(
  "count_matrix_RA.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

dim(counts)
head(counts[, 1:4])


# ----------------------------
# 4. Sample metadata
# ----------------------------

sample_info <- data.frame(
  sample_id = c(
    "SRR4785819",
    "SRR4785820",
    "SRR4785828",
    "SRR4785831",
    "SRR4785979",
    "SRR4785980",
    "SRR4785986",
    "SRR4785988"
  ),
  condition = c(
    "Control",
    "Control",
    "Control",
    "Control",
    "RA",
    "RA",
    "RA",
    "RA"
  )
)

rownames(sample_info) <- sample_info$sample_id
sample_info$condition <- factor(sample_info$condition, levels = c("Control", "RA"))

# Zorg dat de volgorde van de samples in counts overeenkomt met sample_info
counts <- counts[, sample_info$sample_id]

write.csv(sample_info, "results/sample_metadata.csv", row.names = FALSE)


# ----------------------------
# 5. DESeq2 analyse
# ----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = sample_info,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "RA", "Control")
)

summary(res)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_ordered <- res_df[order(res_df$padj), ]

write.csv(
  res_ordered,
  "results/DESeq2_results.csv",
  row.names = TRUE
)


# ----------------------------
# 6. Significante genen
# ----------------------------

sig_genes_padj005 <- subset(
  res_df,
  !is.na(padj) & padj < 0.05
)

write.csv(
  sig_genes_padj005,
  "results/significant_genes_padj005.csv",
  row.names = TRUE
)

# Aantallen up- en downregulated genen
upregulated <- sum(sig_genes_padj005$log2FoldChange > 0)
downregulated <- sum(sig_genes_padj005$log2FoldChange < 0)

upregulated
downregulated


# ----------------------------
# 7. PCA plot
# ----------------------------

vsd <- vst(dds)

pcaData <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 5) +
  theme_bw(base_size = 16) +
  labs(
    title = "PCA Plot: RA vs Control",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

ggsave(
  "figures/PCA_plot.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)


# ----------------------------
# 8. Volcano plot
# ----------------------------

volcano_df <- res_df[!is.na(res_df$padj), ]

volcano_df$group <- "Not significant"

volcano_df$group[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange > 1
] <- "Up regulated"

volcano_df$group[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange < -1
] <- "Down regulated"

top_up <- volcano_df[
  volcano_df$group == "Up regulated",
]

top_up <- top_up[order(top_up$padj), ]
top_up <- head(top_up, 10)

top_down <- volcano_df[
  volcano_df$group == "Down regulated",
]

top_down <- top_down[order(top_down$padj), ]
top_down <- head(top_down, 10)

top_labels <- rbind(top_up, top_down)

p_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = group
  )
) +
  geom_point(size = 1.6, alpha = 0.75) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 3.5,
    max.overlaps = 30
  ) +
  scale_color_manual(
    values = c(
      "Down regulated" = "#4DBBD5",
      "Not significant" = "grey75",
      "Up regulated" = "#E64B35"
    )
  ) +
  theme_bw(base_size = 16) +
  labs(
    title = "Volcano Plot: RA vs Control",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value"),
    color = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

ggsave(
  "figures/Volcano_plot.png",
  p_volcano,
  width = 9,
  height = 7,
  dpi = 300
)


# ----------------------------
# 9. GO analyse
# ----------------------------

gene_symbols <- rownames(sig_genes_padj005)

entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

GO_results <- as.data.frame(ego)

write.csv(
  GO_results,
  "results/GO_results.csv",
  row.names = FALSE
)

p_go <- dotplot(
  ego,
  showCategory = 10,
  label_format = 40,
  font.size = 12
) +
  ggtitle("GO enrichment: Biological Processes") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  "figures/GO_dotplot.png",
  p_go,
  width = 10,
  height = 7,
  dpi = 300
)


# ----------------------------
# 10. KEGG analyse
# ----------------------------

kegg <- enrichKEGG(
  gene = entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

KEGG_results <- as.data.frame(kegg)

write.csv(
  KEGG_results,
  "results/KEGG_results.csv",
  row.names = FALSE
)

p_kegg <- dotplot(
  kegg,
  showCategory = 10,
  label_format = 35,
  font.size = 12
) +
  ggtitle("KEGG pathway enrichment") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  "figures/KEGG_dotplot.png",
  p_kegg,
  width = 10,
  height = 7,
  dpi = 300
)


# ----------------------------
# 11. Pathview - TNF signaling pathway
# ----------------------------
# hsa04668 = TNF signaling pathway

fc_path <- res_df$log2FoldChange
names(fc_path) <- rownames(res_df)

gene_df <- bitr(
  names(fc_path),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

fc_entrez <- fc_path[gene_df$SYMBOL]
names(fc_entrez) <- gene_df$ENTREZID

pathview(
  gene.data = fc_entrez,
  pathway.id = "hsa04668",
  species = "hsa",
  gene.idtype = "entrez",
  out.suffix = "TNF_RA_vs_Control"
)


# ----------------------------
# 12. Session info
# ----------------------------

sink("results/sessionInfo.txt")
sessionInfo()
sink()

message("Analyse afgerond. Controleer de mappen results/ en figures/.")

# Transcriptomics Analyse van Reumatoïde Artritis (RA)

## Projectoverzicht

In dit project is een transcriptomics-analyse uitgevoerd om verschillen in genexpressie tussen patiënten met Reumatoïde Artritis (RA) en gezonde controles te onderzoeken met behulp van RNA-sequencing data.

De analyse omvat:

- Differentiële genexpressieanalyse (DESeq2)
- Principal Component Analysis (PCA)
- Volcano plot visualisatie
- Gene Ontology (GO) enrichment analyse
- KEGG pathway enrichment analyse
- Visualisatie van het TNF-signaleringspad

---

## Dataset

De dataset bestond uit:

- 4 RA-samples
- 4 gezonde controles

De count matrix werd aangeleverd na read alignment en genkwantificatie.

---

## Analyseworkflow

### 1. Datavoorbereiding

- Inlezen van de count matrix
- Opstellen van sample metadata
- Controle van de dataset

### 2. Differentiële genexpressieanalyse

De analyse werd uitgevoerd met het R-package DESeq2.

Selectiecriteria voor significante genen:

- adjusted p-value (padj) < 0,05

De resultaten werden opgeslagen in:

- DESeq2_results.csv
- significant_genes_padj005.csv

### 3. Functionele verrijkingsanalyse

Voor de significante genen werden uitgevoerd:

- GO Biological Process enrichment analyse
- KEGG pathway enrichment analyse

### 4. Visualisaties

De volgende figuren werden gegenereerd:

- PCA Plot
- Volcano Plot
- GO Dotplot
- KEGG Dotplot
- TNF Signaling Pathway

---

## Repositorystructuur

```
Transcriptomics_RA_Final/
│
├── data/
│   ├── count_matrix_RA.csv
│   └── data_raw/
│
├── figures/
│   ├── PCA_plot.png
│   ├── Volcano_plot.png
│   ├── GO_dotplot.png
│   ├── KEGG_dotplot.png
│   └── TNF_signaling_pathway.png
│
├── results/
│   ├── DESeq2_results.csv
│   ├── significant_genes_padj005.csv
│   ├── GO_results.csv
│   ├── KEGG_results.csv
│   └── sample_metadata.csv
│
├── scripts/
│   ├── 01_data_preparation.R
│   └── 02_DESeq2_GO_KEGG_analysis.R
│
├── bronnen/
│
└── README.md
```

---

## Gebruikte software en packages

- R
- DESeq2
- clusterProfiler
- enrichplot
- ggplot2
- pathview

---

## Auteur

Ahmed Abbas
Studentnummer: 5388929

Biologie en Medisch Laboratoriumonderzoek (BML)
NHL Stenden Hogeschool / Van Hall Larenstein
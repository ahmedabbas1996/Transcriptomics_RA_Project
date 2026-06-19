# ============================================================
# 01_FASTQ_mapping_featureCounts.R
# Transcriptomics RA project
#
# Doel:
# Dit script beschrijft de eerste bioinformatica-stappen vanaf FASTQ-bestanden
# tot en met het maken van een count matrix.
#
# Workflow:
#   1. Projectmappen controleren/aanmaken
#   2. Sample metadata maken
#   3. FASTQ-bestanden controleren
#   4. Aantal reads per FASTQ-bestand tellen
#   5. Barplot maken van read counts
#   6. Human reference index bouwen met Rsubread
#   7. Paired-end reads alignen met Rsubread
#   8. Genen tellen met featureCounts
#   9. Count matrix opslaan
#
# Belangrijk:
# De officiële downstream analyse is uitgevoerd met de count_matrix_RA.txt
# die door de docent beschikbaar is gesteld. Dit script laat zien hoe de
# mapping- en counting-stappen technisch zijn uitgevoerd/geoefend.
# ============================================================


# ----------------------------
# 0. Packages
# ----------------------------

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Rsubread", ask = FALSE, update = FALSE)
}

library(Rsubread)


# ----------------------------
# 1. Project folders
# ----------------------------

dir.create("data_raw", showWarnings = FALSE)
dir.create("reference", showWarnings = FALSE)
dir.create("bam", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ----------------------------
# 2. Sample metadata
# ----------------------------

samples <- data.frame(
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

write.csv(
  samples,
  "results/sample_metadata.csv",
  row.names = FALSE
)

samples


# ----------------------------
# 3. FASTQ overview
# ----------------------------

fastq_files <- list.files(
  "data_raw",
  pattern = "\\.fastq$",
  full.names = TRUE
)

fastq_overview <- data.frame(
  file = basename(fastq_files),
  size_MB = round(file.size(fastq_files) / 1024^2, 2)
)

write.csv(
  fastq_overview,
  "results/fastq_overview.csv",
  row.names = FALSE
)

fastq_overview


# ----------------------------
# 4. Count reads in FASTQ files
# ----------------------------
# FASTQ format: 1 read bestaat uit 4 regels.

read_counts <- sapply(
  fastq_files,
  function(f) {
    n_lines <- length(readLines(f))
    n_lines / 4
  }
)

read_summary <- data.frame(
  file = basename(fastq_files),
  reads = as.numeric(read_counts)
)

write.csv(
  read_summary,
  "results/read_counts.csv",
  row.names = FALSE
)

read_summary


# ----------------------------
# 5. Read count barplot
# ----------------------------

png(
  filename = "figures/read_counts_barplot.png",
  width = 1200,
  height = 700
)

barplot(
  read_summary$reads,
  names.arg = read_summary$file,
  las = 2,
  main = "Number of reads per FASTQ file",
  ylab = "Number of reads"
)

dev.off()


# ----------------------------
# 6. Build reference index
# ----------------------------
# Vereiste input:
# reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#
# Let op:
# Deze stap kan lang duren. Sluit RStudio of de laptop niet tijdens buildindex().

reference_file <- "reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
index_basename <- "reference/human_GRCh38_index"

if (!file.exists(reference_file)) {
  stop(
    paste(
      "Reference genome file not found:",
      reference_file,
      "\nPlaats het referentiegenoom in de map reference/."
    )
  )
}

buildindex(
  basename = index_basename,
  reference = reference_file
)


# ----------------------------
# 7. Align paired-end reads
# ----------------------------
# Verwachte FASTQ namen:
# data_raw/SRR4785819_1_subset40k.fastq
# data_raw/SRR4785819_2_subset40k.fastq
# etc.

for (s in samples$sample_id) {

  read_1 <- paste0("data_raw/", s, "_1_subset40k.fastq")
  read_2 <- paste0("data_raw/", s, "_2_subset40k.fastq")
  bam_out <- paste0("bam/", s, ".bam")

  if (!file.exists(read_1)) {
    stop(paste("Missing FASTQ file:", read_1))
  }

  if (!file.exists(read_2)) {
    stop(paste("Missing FASTQ file:", read_2))
  }

  align(
    index = index_basename,
    readfile1 = read_1,
    readfile2 = read_2,
    output_file = bam_out,
    nthreads = 2
  )
}

list.files("bam")


# ----------------------------
# 8. featureCounts
# ----------------------------
# Vereiste input:
# reference/Homo_sapiens.GRCh38.115.gtf.gz

gtf_file <- "reference/Homo_sapiens.GRCh38.115.gtf.gz"

if (!file.exists(gtf_file)) {
  stop(
    paste(
      "GTF annotation file not found:",
      gtf_file,
      "\nPlaats het GTF-bestand in de map reference/."
    )
  )
}

bam_files <- list.files(
  "bam",
  pattern = "\\.bam$",
  full.names = TRUE
)

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  nthreads = 2
)

fc
dim(fc$counts)


# ----------------------------
# 9. Save count matrix
# ----------------------------

counts <- fc$counts

colnames(counts) <- samples$sample_id

write.csv(
  counts,
  "results/gene_counts_from_featureCounts.csv"
)

counts[1:10, ]


# ----------------------------
# 10. Session info
# ----------------------------

sink("results/sessionInfo_mapping_featureCounts.txt")
sessionInfo()
sink()

message("Mapping and featureCounts script finished.")
message("Check results/gene_counts_from_featureCounts.csv")

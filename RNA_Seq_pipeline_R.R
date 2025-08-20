############################################################
# R Script for RNA-seq QC and Preprocessing
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR","limma","tximport"), ask = FALSE, update = FALSE)

# Step 1: Quality Control (FastQC) and Trimming (fastp)
############################################################

# Load libraries for downstream (only needed later, not for system calls now)
if(!require(edgeR)) install.packages("edgeR", repos="http://cran.us.r-project.org")
if(!require(tximport)) BiocManager::install("tximport")
library(edgeR)
library(limma)
library(tximport)
# Define paths
raw_dir <- "raw_data"          # Folder with raw fastq.gz files
qc_dir <- "QC_reports"         # Folder for FastQC output
trim_dir <- "trimmed_data"     # Folder for trimmed reads

# Create directories if they don't exist
dir.create(qc_dir, showWarnings = FALSE)
dir.create(trim_dir, showWarnings = FALSE)

# List raw FASTQ files
fastq_files <- list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Step 1a: Run FASTQC on raw reads
for(fq in fastq_files){
  system(paste("fastqc -o", qc_dir, fq))
}

# Step 1b: Run fastp for trimming (default parameters)
# Output paired trimmed reads into trim_dir
for(i in seq(1, length(fastq_files), by=2)){   # assuming paired-end reads
  fq1 <- fastq_files[i]
  fq2 <- fastq_files[i+1]
  sample_name <- sub("_R1.*", "", basename(fq1))
  
  out1 <- file.path(trim_dir, paste0(sample_name, "_R1_trimmed.fastq.gz"))
  out2 <- file.path(trim_dir, paste0(sample_name, "_R2_trimmed.fastq.gz"))
  html <- file.path(trim_dir, paste0(sample_name, "_fastp.html"))
  json <- file.path(trim_dir, paste0(sample_name, "_fastp.json"))
  
  cmd <- paste("fastp -i", fq1, "-I", fq2,
               "-o", out1, "-O", out2,
               "-h", html, "-j", json)
  system(cmd)
}

# Step 1c: Run FASTQC again on trimmed reads
trimmed_files <- list.files(trim_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
for(fq in trimmed_files){
  system(paste("fastqc -o", qc_dir, fq))
}

############################################################
# Output: 
# - QC reports (before & after trimming) in QC_reports/
# - Cleaned FASTQ files in trimmed_data/
############################################################
############################################################
# R Script for De novo Transcriptome Assembly
# Step 2: Trinity assembly + assembly-stats evaluation
############################################################

# Define directories
trim_dir <- "trimmed_data"       # from Step 1
assembly_dir <- "trinity_assembly"
stats_dir <- "assembly_stats"

# Create directories if they don't exist
dir.create(assembly_dir, showWarnings = FALSE)
dir.create(stats_dir, showWarnings = FALSE)

# List trimmed FASTQ files (paired-end)
trimmed_files <- list.files(trim_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Make sure files are sorted properly (R1 before R2)
fq1 <- trimmed_files[grep("_R1", trimmed_files)]
fq2 <- trimmed_files[grep("_R2", trimmed_files)]

# Run Trinity (default parameters)
# Adjust --max_memory and --CPU based on your system
cmd_trinity <- paste(
  "Trinity",
  "--seqType fq",
  "--left", paste(fq1, collapse=","),
  "--right", paste(fq2, collapse=","),
  "--CPU 8",                        # change based on available cores
  "--max_memory 50G",               # change based on your RAM
  "--output", assembly_dir
)
system(cmd_trinity)

# The assembled transcriptome will be saved as:
# trinity_assembly/Trinity.fasta

# Step 2b: Assembly quality evaluation with assembly-stats
trinity_fasta <- file.path(assembly_dir, "Trinity.fasta")
stats_out <- file.path(stats_dir, "stats.txt")

cmd_stats <- paste("assembly-stats", trinity_fasta, ">", stats_out)
system(cmd_stats)

############################################################
# Output:
# - Assembled transcriptome: trinity_assembly/Trinity.fasta
# - Assembly statistics: assembly_stats/stats.txt
############################################################
############################################################
# R Script for Transcript Quantification and DEG Analysis
# Step 3: RSEM quantification + edgeR DEG
############################################################

# Load Bioconductor libraries
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("edgeR", "limma"), ask=FALSE, update=FALSE)

library(edgeR)
library(limma)

# Paths
assembly_dir <- "trinity_assembly"
trim_dir <- "trimmed_data"
rsem_dir <- "rsem_results"
deg_dir <- "DEG_results"

dir.create(rsem_dir, showWarnings = FALSE)
dir.create(deg_dir, showWarnings = FALSE)

# Step 3a: Prepare reference for RSEM (only once)
trinity_fasta <- file.path(assembly_dir, "Trinity.fasta")
rsem_ref <- file.path(rsem_dir, "Trinity_RSEM")

cmd_prepare <- paste("rsem-prepare-reference",
                     "--transcript-to-gene-map", "transcript_to_gene.map", # optional if available
                     trinity_fasta, rsem_ref)
system(cmd_prepare)

# Step 3b: Quantify expression for each sample using RSEM
# Loop over pairs of trimmed fastq files
trimmed_files <- list.files(trim_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

fq1 <- trimmed_files[grep("_R1", trimmed_files)]
fq2 <- trimmed_files[grep("_R2", trimmed_files)]

for(i in seq_along(fq1)){
  sample_name <- sub("_R1.*", "", basename(fq1[i]))
  out_prefix <- file.path(rsem_dir, sample_name)
  
  cmd_rsem <- paste("rsem-calculate-expression",
                    "--paired-end",
                    "-p 8",                        # adjust CPU
                    fq1[i], fq2[i],
                    rsem_ref, out_prefix)
  system(cmd_rsem)
}

# At this point, RSEM will produce gene and isoform expression results:
# *.genes.results and *.isoforms.results for each sample

# Step 3c: Import counts into R for edgeR
gene_files <- list.files(rsem_dir, pattern="genes.results$", full.names=TRUE)
sample_names <- sub("\\.genes.results", "", basename(gene_files))

# Read counts into DGEList
counts_list <- lapply(gene_files, function(f) read.delim(f, row.names=1))
count_matrix <- do.call(cbind, lapply(counts_list, function(x) x$expected_count))
colnames(count_matrix) <- sample_names

dge <- DGEList(counts=count_matrix)
dge <- calcNormFactors(dge)

# Step 3d: Experimental groups
# Define your groups (edit names to match actual sample names)
group <- factor(c("SA", "SA-C", "SA-P", "SA-E"))  
dge$samples$group <- group

# Step 3e: Filter out uniquely expressed genes (only expressed in one sample)
keep <- rowSums(cpm(dge) > 1) > 1   # keep genes expressed in >1 sample
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Step 3f: Create design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Step 3g: Estimate dispersions
dge <- estimateDisp(dge, design)

# Step 3h: Fit GLM and test
fit <- glmFit(dge, design)

# Comparisons: SA as control
contrasts <- makeContrasts(
  SA_CvsSA = `SA-C` - SA,
  SA_PvsSA = `SA-P` - SA,
  SA_EvsSA = `SA-E` - SA,
  levels=design
)

# Run LRT for each comparison
lrt_SA_C <- glmLRT(fit, contrast=contrasts[,"SA_CvsSA"])
lrt_SA_P <- glmLRT(fit, contrast=contrasts[,"SA_PvsSA"])
lrt_SA_E <- glmLRT(fit, contrast=contrasts[,"SA_EvsSA"])

# Step 3i: Save DEG results
write.table(topTags(lrt_SA_C, n=Inf), file=file.path(deg_dir, "DEG_SA-C_vs_SA.txt"), sep="\t", quote=FALSE)
write.table(topTags(lrt_SA_P, n=Inf), file=file.path(deg_dir, "DEG_SA-P_vs_SA.txt"), sep="\t", quote=FALSE)
write.table(topTags(lrt_SA_E, n=Inf), file=file.path(deg_dir, "DEG_SA-E_vs_SA.txt"), sep="\t", quote=FALSE)

############################################################
# Output:
# - RSEM expression tables: rsem_results/
# - DEG results: DEG_results/ (per comparison)
# - Uniquely expressed genes excluded automatically
############################################################
############################################################
# R Script for Functional Annotation of DEGs
# Step 4: TransDecoder + BLAST + annotation filtering
############################################################

# Define directories
assembly_dir <- "trinity_assembly"
annotation_dir <- "annotation_results"
blast_dir <- file.path(annotation_dir, "BLAST")

dir.create(annotation_dir, showWarnings = FALSE)
dir.create(blast_dir, showWarnings = FALSE)

# Step 4a: Run TransDecoder to extract coding regions
trinity_fasta <- file.path(assembly_dir, "Trinity.fasta")

cmd_transdecoder_longorfs <- paste("TransDecoder.LongOrfs -t", trinity_fasta)
cmd_transdecoder_predict  <- paste("TransDecoder.Predict -t", trinity_fasta)

system(cmd_transdecoder_longorfs)
system(cmd_transdecoder_predict)

# Output:
# - Trinity.fasta.transdecoder_dir/
# - Trinity.fasta.transdecoder.pep  (predicted proteins)

# Step 4b: Run BLASTp against UniProt (Swiss-Prot reviewed)
pep_file <- paste0(trinity_fasta, ".transdecoder.pep")
blast_out <- file.path(blast_dir, "blastp_results.txt")

# NOTE: uniprot_sprot.fasta should be downloaded locally and formatted with makeblastdb
# Example: makeblastdb -in uniprot_sprot.fasta -dbtype prot

cmd_blast <- paste("blastp -query", pep_file,
                   "-db uniprot_sprot.fasta",
                   "-out", blast_out,
                   "-evalue 1e-5 -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen evalue bitscore'")
system(cmd_blast)

# Step 4c: Import BLAST results and filter by ≥50% similarity
blast_results <- read.delim(blast_out, header=FALSE, sep="\t")
colnames(blast_results) <- c("QueryID","SubjectID","PercentIdentity",
                             "AlignmentLength","Mismatches","GapOpens",
                             "Evalue","Bitscore")

# Apply ≥50% identity filter
filtered_hits <- subset(blast_results, PercentIdentity >= 50)

# Save filtered results
write.table(filtered_hits, file=file.path(blast_dir, "blastp_filtered_50.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

############################################################
# Output:
# - Predicted proteins: Trinity.fasta.transdecoder.pep
# - Raw BLAST results: BLAST/blastp_results.txt
# - Filtered hits (≥50% identity): BLAST/blastp_filtered_50.txt
############################################################


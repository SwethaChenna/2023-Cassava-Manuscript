library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)
library(collections)
library(TranscriptomeReconstructoR)
library(devtools)

devtools::source_url("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/batch_read_track_data.R?raw=TRUE")

# Prepare Seqinfo object ------------------------------------------------------

si <- file.path("path_to_file", "infile.txt") %>% read_tsv(col_names = c("chr", "len"), col_types = "ci") ### file with chromosome names and lengths in two columns
si <- Seqinfo(seqnames = si$chr, seqlengths = si$len)

##### Prepare input data -----------------------------------------------------

# ONT data from Bam files:
ONT_wt <- load_BAM_files("path_to_bam_file.bam", mode = "long_read")
ONT_cold <- load_BAM_files("path_to_bam_file.bam", mode = "long_read")
saveRDS(ONT_wt, "ONT_wt.RDS")
saveRDS(ONT_cold, "ONT_Cold.RDS")

# STRIPE-Seq from Bam files:
tss_wt_files <- c("path_to_bam_file_rep1.bam", 
                  "path_to_bam_file_rep2.bam", 
                  "path_to_bam_file_rep3.bam", 
                  "path_to_bam_file_rep4.bam")
tss_cold_files <- c("path_to_bam_file_rep1.bam",
                    "path_to_bam_file_rep2.bam", 
                    "path_to_bam_file_rep3.bam", 
                    "path_to_bam_file_rep4.bam")
tss_wt <- load_BAM_files(tss_wt_files, mode = "tss")
tss_cold <- load_BAM_files(tss_cold_files, mode = "tss")
saveRDS(tss_wt, "tss_wt.RDS")
saveRDS(tss_cold, "tss_cold.RDS")

# Quant-Seq from Bam files:
pas_wt_files <- c("path_to_bam_file_rep1.bam",
                  "path_to_bam_file_rep2.bam", 
                  "path_to_bam_file_rep3.bam", 
                  "path_to_bam_file_rep4.bam")
pas_cold_files <- c("path_to_bam_file_rep1.bam",
                    "path_to_bam_file_rep2.bam", 
                    "path_to_bam_file_rep3.bam", 
                    "path_to_bam_file_rep4.bam")
pas_wt <- load_BAM_files(pas_wt_files, mode = "pas")
pas_cold <- load_BAM_files(pas_cold_files, mode = "pas")
saveRDS(pas_wt, "pas_wt.RDS")
saveRDS(pas_cold, "pas_cold.RDS")

# ssRNA-Seq data (as a substitute for plaNET-seq) from Bam files:
ssrna_wt_files <- c("path_to_bam_file_rep1.bam",
                    "path_to_bam_file_rep2.bam", 
                    "path_to_bam_file_rep3.bam", 
                    "path_to_bam_file_rep4.bam")
ssrna_cold_files <- c("path_to_bam_file_rep1.bam",
                      "path_to_bam_file_rep2.bam", 
                      "path_to_bam_file_rep3.bam", 
                      "path_to_bam_file_rep4.bam")
ssrna_wt <- load_BAM_files(ssrna_wt_files, ngs_mode = "PE")
ssrna_cold <- load_BAM_files(ssrna_cold_files, ngs_mode = "PE")
saveRDS(ssrna_wt, "ssrna_wt.RDS")
saveRDS(ssrna_cold, "ssrna_cold.RDS")


##### Load the prepared data ----------------------------------------

ont_wt_org <- readRDS("ONT_wt.RDS")
ont_cold <- readRDS("ONT_Cold.RDS")
tss_data <- readRDS("tss_data.RDS")
tss_cold <- tss_data[1:4]
tss_wt <- tss_data[5:8]
pas_data <- readRDS("pas_data.RDS")
pas_cold <- pas_data[1:4]
pas_wt <- pas_data[5:8]
rna_wt <- readRDS("ssrna_wt_data.RDS")
#rna_cold <- readRDS("ssrna_cold_data.RDS")


##### Call TSS and PAS ----------------------------------------------

tss_tc_wt <- call_TCs(tss_wt, min_support = 3, min_tpm = 0.5)
tss_tc_cold <- call_TCs(tss_cold, min_support = 3, min_tpm = 0.5)
pas_tc_wt <- call_TCs(pas_wt, min_support = 3, min_tpm = 0.1)
#pas_tc_cold <- call_TCs(pas_cold, min_support = 3, min_tpm = 0.1)

export(tss_tc_wt, "TSS_TC_wt.bed", format = "BED")
export(tss_tc_cold, "TSS_TC_cold.bed", format = "BED")
export(pas_tc_wt, "PAS_TC_wt.bed", format = "BED")
export(pas_tc_cold, "PAS_TC_cold.bed", format = "BED")

saveRDS(tss_tc_wt, "tss_tc_wt.RDS")
saveRDS(tss_tc_cold, "tss_tc_cold.RDS")
saveRDS(pas_tc_wt, "pas_tc_wt.RDS")
saveRDS(pas_tc_cold, "pas_tc_cold.RDS")

tss_tc_wt <- readRDS("tss_tc_wt.RDS")
tss_tc_cold <- readRDS("tss_tc_cold.RDS")
pas_tc_wt <- readRDS("pas_tc_wt.RDS")
pas_tc_cold <- readRDS("pas_tc_cold.RDS")


##### Run TrRR on wt data ---------------------------------------------

# Correct long reads:
ont_wt <- extend_long_reads_to_TSS_and_PAS(ont_wt, tss_tc_wt, pas_tc_wt)
ont_wt <- adjust_exons_of_long_reads(ont_wt)
ont_wt <- detect_alignment_errors(ont_wt)

# Call de novo gene and transcript model:
out_wt <- call_transcripts_and_genes(ont_wt)
hml_genes_wt <- out_wt[[1]]
hml_tx_wt <- out_wt[[2]]
fusion_genes_wt <- out_wt[[3]]
fusion_tx_wt <- out_wt[[4]]
reads_free_wt <- out_wt[[5]]

# Also call transcription units from RNA-seq:
trans_wt <- call_transcribed_intervals(rna_wt)
transcribed_wt <- trans_wt[[1]]
gaps_wt <- trans_wt[[2]]

export(transcribed_wt, "Transcribed_intervals_RNAseq_wt.bed", format = "BED")

results_wt <- process_nascent_intervals(hml_genes_wt, transcribed_wt, tss_tc_wt, pas_tc_wt, reads_free_wt, gaps_wt)
hml_genes_v2_wt <- results_wt[[1]]
lowexpr_wt <- results_wt[[3]]

rtracklayer::export(hml_genes_v2_wt, "Called_genes_wt.bed", format = "BED")
write_grl_as_bed12(hml_tx_wt, "Called_transcripts_wt.bed")
rtracklayer::export(fusion_genes_wt, "Fusion_genes_wt.bed", format = "BED")
write_grl_as_bed12(fusion_tx_wt, "Fusion_transcripts_wt.bed")
rtracklayer::export(lowexpr_wt, "Low_expressed_genes_wt.bed", format = "BED")

# Also save all results as RDS file:
final_out_wt <- list("hml_genes_wt" = hml_genes_wt, "hml_genes_v2_wt" = hml_genes_v2_wt, "hml_tx_wt" = hml_tx_wt, "fusion_genes_wt" = fusion_genes_wt, "fusion_tx_wt" = fusion_tx_wt, "lowexp_wt" = lowexpr_wt)
saveRDS(final_out_wt, "final_out_wt.RDS")


##### Also run TrRR on cold_3h data -----------------------------------------

# Correct long reads:
ont_cold <- extend_long_reads_to_TSS_and_PAS(ont_cold, tss_tc_cold, pas_tc_cold)
ont_cold <- adjust_exons_of_long_reads(ont_cold)
ont_cold <- detect_alignment_errors(ont_cold)

# Call de novo gene and transcript model:
out_cold <- call_transcripts_and_genes(ont_cold)
hml_genes_cold <- out_cold[[1]]
hml_tx_cold <- out_cold[[2]]
fusion_genes_cold <- out_cold[[3]]
fusion_tx_cold <- out_cold[[4]]
reads_free_cold <- out_cold[[5]]

# Also call transcription units from RNA-seq:
trans_cold <- call_transcribed_intervals(rna_cold)
transcribed_cold <- trans_cold[[1]]
gaps_cold <- trans_cold[[2]]

export(transcribed_cold, "Transcribed_intervals_RNAseq_cold.bed", format = "BED")

results_cold <- process_nascent_intervals(hml_genes_cold, transcribed_cold, tss_tc_cold, pas_tc_cold, reads_free_cold, gaps_cold)
hml_genes_v2_cold <- results_cold[[1]]
lowexpr_cold <- results_cold[[3]]

rtracklayer::export(hml_genes_v2_cold, "Called_genes_cold3h.bed", format = "BED")
write_grl_as_bed12(hml_tx_cold, "Called_transcripts_cold3h.bed")
rtracklayer::export(fusion_genes_cold, "Fusion_genes_cold3h.bed", format = "BED")
write_grl_as_bed12(fusion_tx_cold, "Fusion_transcripts_cold3h.bed")
rtracklayer::export(lowexpr_cold, "Low_expressed_genes_cold3h.bed", format = "BED")

# Also save all results as RDS file:
final_out_cold <- list("hml_genes_cold" = hml_genes_cold, "hml_genes_v2_cold" = hml_genes_v2_cold, "hml_tx_cold" = hml_tx_cold, "fusion_genes_cold" = fusion_genes_cold, "fusion_tx_cold" = fusion_tx_cold, "lowexp_cold" = lowexpr_cold)
saveRDS(final_out_cold, "final_out_cold.RDS")


library(GEOquery)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(GenomicAlignments)

# This was the code used to generate all gviz plots used in the presentation and report



# ----------------- Getting hnRNPM binding CLIP data- --------------------------

chromosome <- "chr3"       # Set chromosome number

bed_file_path <- "/mnt/gtklab01/darren/GSE227047/GSE227048_HMLE_hnM_sig_peaks.bed.gz"
hnRNPM_peaks <- import(bed_file_path, format = "BED")

# Subset for a specific chromosome
hnRNPM_peaks_ch <- hnRNPM_peaks[ chromosome == seqnames(hnRNPM_peaks)]
# Create annotation track for hnRNPM peaks
hnRNPM_track <- AnnotationTrack(hnRNPM_peaks_ch, name = "hnRNPM Peaks", genome = "hg38", stackHeight =1)


# ------------------------- Getting RNAseq Tracks ------------------------------

ylim <- c(0,1000)   # height range to visualise RNAseq peaks

# Get RNAseq data
wt_rna_seq_path1 <- "/mnt/gtklab01/darren/GSE227047_REDO/STAR_CryEx_trimmed/GSM7090417_Aligned.sortedByCoord.out.bam"
wt_rna_seq_path2 <- "/mnt/gtklab01/darren/GSE227047_REDO/STAR_CryEx_trimmed/GSM7090419_Aligned.sortedByCoord.out.bam"
wt_rna_seq_path3 <- "/mnt/gtklab01/darren/GSE227047_REDO/STAR_CryEx_trimmed/GSM7090421_Aligned.sortedByCoord.out.bam"

ko_rna_seq_path1 <- "/mnt/gtklab01/darren/GSE227047_REDO/STAR_CryEx_trimmed/GSM7090418_Aligned.sortedByCoord.out.bam"
ko_rna_seq_path2 <- "/mnt/gtklab01/darren/GSE227047_REDO/STAR_CryEx_trimmed/GSM7090420_Aligned.sortedByCoord.out.bam"
ko_rna_seq_path3 <- "/mnt/gtklab01/darren/GSE227047_REDO/STAR_CryEx_trimmed/GSM7090422_Aligned.sortedByCoord.out.bam"

wt_align_track_wt1 <- AlignmentsTrack(wt_rna_seq_path1, isPaired = TRUE, genome = "hg38", name = "WT RNA-seq1", ylim = ylim)
wt_align_track_wt2 <- AlignmentsTrack(wt_rna_seq_path2, isPaired = TRUE, genome = "hg38", name = "WT RNA-seq2", ylim = ylim)
wt_align_track_wt3 <- AlignmentsTrack(wt_rna_seq_path3, isPaired = TRUE, genome = "hg38", name = "WT RNA-seq3", ylim = ylim)

ko_align_track_ko1 <- AlignmentsTrack(ko_rna_seq_path1, isPaired = TRUE, genome = "hg38", name = "KO RNA-seq1", ylim = ylim)
ko_align_track_ko2 <- AlignmentsTrack(ko_rna_seq_path2, isPaired = TRUE, genome = "hg38", name = "KO RNA-seq2", ylim = ylim)
ko_align_track_ko3 <- AlignmentsTrack(ko_rna_seq_path2, isPaired = TRUE, genome = "hg38", name = "KO RNA-seq3", ylim = ylim)

# Overlay the RNAseq tracks and create the annotation tracks
wt_align_track_overlay <- OverlayTrack(trackList=list(wt_align_track_wt1,wt_align_track_wt2,wt_align_track_wt3), name = "WT RNA-seq")
ko_align_track_overlay <- OverlayTrack(trackList=list(ko_align_track_ko2,ko_align_track_ko2,ko_align_track_ko3), name = "KO RNA-seq")

# ----------------------- Getting Genome Axis Track ----------------------------

genome_axis <- GenomeAxisTrack()

# ------------- Creating Reference and Cryptic Gene Model Tracks ---------------

# getting the Reference genome assembly .gtf file
txdb <- makeTxDbFromGFF("/mnt/gtklab01/linglab/hsapiens_annotation_files/gencode.v40.primary_assembly.annotation.gtf", format="gtf")
# getting the genome assembly merged with cryptic exons file
txdb_cryptic_exon <- makeTxDbFromGFF("/mnt/cbis/home/e0958289/GSE227047/SpliCeAT_final/augment_transcriptome/results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf", format="gtf")

# Creating the gene model tracks
gene_track <- GeneRegionTrack(txdb, genome = "hg38", name = "Transcripts", chromosome = chromosome)
gene_track_cryptic <- GeneRegionTrack(txdb_cryptic_exon, genome = "hg38", name = "Cryptic exons", chromosome = chromosome, fill = "red")

# creating a sequence track (if needed)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#library(BSgenome.Hsapiens.UCSC.hg38)
#sTrack <- SequenceTrack(Hsapiens)


# ------------------------------ Plotting Tracks -------------------------------

options(ucscChromosomeNames = FALSE)

start <- 100519533
end <- 100519608


ht <- HighlightTrack(trackList = list(genome_axis, hnRNPM_track, wt_align_track_overlay, ko_align_track_overlay, gene_track, gene_track_cryptic),
                     #start = 35610620,
                     #end = 35640368,
                     chromosome = chromosome)


# ------------------- Top 3 Up regulatd genes -------------------


# ITGB1 integrin subunit beta 1, cryptc exon at chr10 32915810 - 32915842
# ylim 0-25000

plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.075, 0.075, 0.525, 0.525, 0.2, 0.07),
           from = start - 4800 ,
           to = end + 5000)

# "TYW3 tRNA-yW synthesizing protein 3 homolog Retained intron at chr1 74744972 74746864
# ylim 0-1700, there are no hnRNPM binding sites

plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.00, 0.45, 0.45, 0.2, 0.1),
           from = start - 10000 ,
           to = end + 8000)


# TRAPPC10 trafficking protein particle complex subunit 10 chr21 
# retained intron at 44038176 44052279
# ylim 0-1000, there are no hnRNPM binding sites

plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.5, 0.5, 0.2, 0.15),
           from = start - 4000 ,
           to = end + 4000)



# ------------------- Top 3 Down regulatd genes -------------------


# CABIN1 calcineurin binding protein 1 Cryptic exon at chr22:24068467-24070799
# ylim 0-1000, there are no hnRNPM binding sites

plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.5, 0.5, 0.2,0.5),
           from = start - 3000 ,
           to = end + 6000)


# MTMR3 myotubularin related protein 3, overall reduced expression, chr22 30003898-30004021
# ylim 0-1000, there are no hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.5, 0.5, 0.2, 0.1),
           from = start - 3000 ,
           to = end + 6000)


# DCUN1D4 defective in cullin neddylation 1 domain containing 4 chr4:51908976-51908979
# Lower expression
# ylim 0-1000, there are no hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.5, 0.5, 0.2, 0.2),
           from = start - 6000 ,
           to = end + 6000)



# --------------- Replicating the plots published by the paper -----------------

# ACYP1 acylphosphatase 1 chr14:75056447-75056560
# ylim 0-1500, there are no hnRNPM binding sites

plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 50,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.00, 0.5, 0.5, 0.2, 0.1),
           from = start - 4500 ,
           to = end + 8000)


# TMEM45A transmembrane protein 45A, erroneous start site, chr3 100519533-100519608
# ylim 0-1000 there are no hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           #coverageHeight = 2,
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.5, 0.5, 0.2,0.1),
           from = start - 500 ,
           to = end + 300)


# RBM34 RNA binding motif protein 34 chr1 235143711 235145015
# ylim 0-1500, there are no hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.5, 0.5, 0.2, 0.2),
           from = start - 7000 ,
           to = end + 5000)


# ----------------- Genes missed by CryEx but detected by SpliCeAT ----------------

# OFD1 centriole and centriolar satellite protein chrX:13739236-13739314
# ylim 0-1000, almost no difference, some hnRNPM binding sites but no Cryptic exon
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.1, 0.5, 0.5, 0.2, 0.1),
           from = start - 3000 ,
           to = end + 8000)


# FBH1 F-box DNA helicase 1 chr10:5891019-5891116
# FBH1 not picked up by CryEx, very minor differential splicing w low RNAseq levels
# ylim 0-1000, there are no hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.55, 0.55, 0.15, 0.15),
           from = start - 1000 ,
           to = end + 10000)


# CRYL1 crystallin lambda 1 chr13:20426695-20427305
# Region highlight was not flagged by CryEx, but flagged by Spliceat
# was high confidence but overall no significant differential splicing 
# ylim 0,100, there are no hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.0, 0.55, 0.55, 0.15, 0.15),
           from = start - 2000 ,
           to = end + 6000)


# ----------------- Cryptic exons found by CryEx but not SpliCeAT ----------------

# IMPAD1	chr8:57880494-57880557	
# NOTHING
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.1, 0.5, 0.5, 0.2, 0.15),
           from = start - 10000 ,
           to = end + 10000)


# GPR149	chr3 154079069-154079129	154078385	154082516
# NOTHING
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.1, 0.5, 0.5, 0.2, 0.15),
           from = start - 10000 ,
           to = end + 10000)



# MRE11A	chr11:94194404-94194595	94194202	94194945
# MRE11A	chr11:94195046-94195233	94194202	94197278
# NOTHING except for some hnRNPM binding sites
plotTracks(list(ht),
           transcriptAnnotation = "symbol",
           type = c("coverage", "sashimi"),
           #type = c("coverage"),
           #coverageHeight = 2,
           sashimiScore = 20,
           sashimiNumbers=TRUE,
           sizes = c(0.1, 0.1, 0.5, 0.5, 0.2, 0.15),
           from = start - 2000 ,
           to = end + 3000)


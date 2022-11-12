library(tidyverse)
library(ChIPseeker)
library(rtracklayer)
library(optparse)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
# setwd("script")
# option_list <- list(
#     make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input file")
#     )
# opt <- parse_args(OptionParser(option_list = option_list))
# dir_name <- basename(opt$input) %>% str_extract(".*(?=\\.ucsc.*)") %>% paste0("../data/06_annotation/", .)
# create ../data/06_annotation/ if not exist
# if (!dir.exists("./data/06_annotation") dir.create("./data/06_annotation")
# To import broadPeak files
extraCols_broadPeak <- c(
    signalValue = "numeric", pValue = "numeric",
    qValue = "numeric"
)
gr_broadPeak <- import("../data/05_peak_calling/H3K4me1_B.ucsc.bam_peaks.broadPeak",
    format = "BED",
    extraCols = extraCols_broadPeak
)
# Overall annotation status
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
anno <- annotatePeak(gr_broadPeak, tssRegion = c(-3000, 3000), TxDb = txdb,annoDb = "org.Mm.eg.db")
# Visulization
covplot(gr_broadPeak, ylab = "Coverage", xlab = "Genomic Position")
peakHeatmap(gr_broadPeak, TxDb = txdb, upstream = 3000, downstream = 3000, title = "H3K4me1_B")
plotAvgProf2(gr_broadPeak,
    TxDb = txdb, upstream = 3000, downstream = 3000,
    ylab = "Read count frequency", xlab = "Genomic Region (5'->3')",
    conf = 0.95, resample = 1000
)
plotAnnoPie(anno)
plotAnnoBar(anno)
upsetplot(anno,vennpie = TRUE)

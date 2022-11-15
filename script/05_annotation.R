library(tidyverse)
library(ChIPseeker)
library(rtracklayer)
library(optparse)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(parallel)
library(foreach)
library(doParallel)
# add option mc.cores = 40L
options(mc.cores = 40L)
# setwd("script")
# option_list <- list(
#     make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input file")
#     )
# opt <- parse_args(OptionParser(option_list = option_list))
# dir_name <- basename(opt$input) %>% str_extract(".*(?=\\.ucsc.*)") %>% paste0("../data/06_annotation/", .)
# create ../data/06_annotation/ if not exist
if (!dir.exists("../data/06_annotation")) dir.create("../data/06_annotation")
# To import broadPeak files
extraCols_broadPeak <- c(
    signalValue = "numeric", pValue = "numeric",
    qValue = "numeric"
)
gr_list <- list()
for (file in list.files("../data/05_peak_calling", pattern = "broadPeak", full.names = TRUE)) {
    gr <- import(file, format = "BED", extraCols = extraCols_broadPeak)
    gr_list[[basename(file) %>% str_extract(".*(?=\\.ucsc.*)")]] <- gr
}
# Overall annotation status
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
anno_list <- mclapply(gr_list, annotatePeak,
    tssRegion = c(-3000, 3000), TxDb = txdb
)
# Visulization
# use for loop to traverse every four elements in gr_list, and use covplot to plot
myCluster <- makeCluster(50)
registerDoParallel(myCluster)
foreach(i = seq(1, length(gr_list), by = 4), .inorder = TRUE, .packages = c("ChIPseeker","tidyverse")) %dopar% {
    pdf(paste0("../data/06_annotation/", names(gr_list)[i], "_cov",".pdf"))
    covplot(gr_list[i:(i + 3)], ylab = "Coverage", xlab = "Genomic Position") + facet_grid(chr ~ .id)
    dev.off()
}
for (i in seq(1, length(gr_list), by = 4)) {
    pdf(paste0("../data/06_annotation/", names(gr_list)[i], "_heatmap", ".pdf"))
    peakHeatmap(gr_list[i:(i + 3)], color = rainbow(n = 4), TxDb = txdb, upstream = 3000, downstream = 3000) +
        facet_grid(chr ~ .id)
    dev.off()
}
foreach(i = seq(1, length(gr_list), by = 4), .inorder = TRUE, .packages = c("ChIPseeker","tidyverse")) %dopar% {
    pdf(paste0("../data/06_annotation/", names(gr_list)[i], "_avgprof", ".pdf"))
    plotAvgProf2(gr_list[i:(i + 3)],
        TxDb = txdb, upstream = 3000, downstream = 3000, ylab = "Read count frequency", xlab = "Genomic Region (5'->3')",
        conf = 0.95, resample = 1000
    ) + facet_grid(chr ~ .id)
    dev.off()
}
foreach(i = seq(1, length(anno_list), by = 1), .inorder = TRUE, .packages = c("ChIPseeker","tidyverse")) %dopar% {
    pdf(paste0("../data/06_annotation/", names(anno_list)[i], "_annopie", ".pdf"))
    plotAnnoPie(anno_list[[i]])
    dev.off()
}
foreach(i = seq(1, length(anno_list), by = 1), .inorder = TRUE, .packages = c("ChIPseeker","tidyverse")) %dopar% {
    pdf(paste0("../data/06_annotation/", names(anno_list)[i], "_annobar", ".pdf"))
    plotAnnoBar(anno_list[i])
    dev.off()
}
foreach(i = seq(1, length(anno_list), by = 1), .inorder = TRUE, .packages = c("ChIPseeker","tidyverse")) %dopar% {
    pdf(paste0("../data/06_annotation/", names(anno_list)[i], "upset", ".pdf"))
    upsetplot(anno_list[[i]], vennpie = TRUE)
    dev.off()
}
stopImplicitCluster()

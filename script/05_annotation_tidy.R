library(tidyverse)
library(ChIPseeker)
library(ggpubr)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(furrr)
library(progressr)
handlers(global = TRUE)
options(future.globals.maxSize = 891289600)
# create ../data/06_annotation/ if not exist
if (!dir.exists("../data/06_annotation")) dir.create("../data/06_annotation")
if (!dir.exists("../data/06_annotation_tidy")) dir.create("../data/06_annotation_tidy")
if (!dir.exists("../data/06_annotation_tidy/annobar")) dir.create("../data/06_annotation_tidy/annobar")
if (!dir.exists("../data/06_annotation_tidy/annopie")) dir.create("../data/06_annotation_tidy/annopie")
if (!dir.exists("../data/06_annotation_tidy/avgprof")) dir.create("../data/06_annotation_tidy/avgprof")
if (!dir.exists("../data/06_annotation_tidy/cov")) dir.create("../data/06_annotation_tidy/cov")
if (!dir.exists("../data/06_annotation_tidy/heatmap")) dir.create("../data/06_annotation_tidy/heatmap")
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
# get all cell types in gr_list names
# cell_types <- names(gr_list) %>%
#     str_extract("(?<=[:alnum:]_).*") %>%
#     unique() %>%
#     factor() %>%
#     relevel(c(
#         "LT_HSC", "ST_HSC", "MPP", "CMP", "GMP", "MF",
#         "GN", "Mono", "CLP", "B", "CD4", "CD8", "NK", "MEP", "EryA", "EryB"
#     ))
cell_types <- c(
    "LT_HSC", "ST_HSC", "MPP", "CMP", "GMP", "MF",
    "GN", "Mono", "CLP", "B", "CD4", "CD8", "NK", "MEP", "EryA", "EryB"
)
# get all histone modifications in gr_list names
histone_mod <- names(gr_list) %>%
    str_extract("[:alnum:]*(?=_.*)") %>%
    unique()
# Overall annotation status
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
anno_list <- lapply(gr_list, function(x) annotatePeak(x, tssRegion = c(-3000, 3000), TxDb = txdb))
package_export <- c("tidyverse", "ChIPseeker")
# Visulization
# use for loop to traverse every four elements in gr_list, and use covplot to plot
plan(multisession, workers = 4)
foreach(i = seq(1, length(gr_list), by = 4), .packages = package_export) %dopar% {
    tmp <- covplot(gr_list[i:(i + 3)], ylab = "Coverage", xlab = "Genomic Position") + facet_grid(chr ~ .id)
    ggsave(tmp, filename = paste0("../data/06_annotation/", names(gr_list)[i], "_cov", ".pdf"), width = 15, height = 10)
}
foreach(i = seq(1, length(gr_list), by = 4), .packages = package_export, .export = txdb) %dopar% {
    pdf(paste0("../data/06_annotation/", names(gr_list)[i], "_heatmap", ".pdf"))
    peakHeatmap(gr_list[i:(i + 3)], color = rainbow(n = 4), TxDb = txdb, upstream = 3000, downstream = 3000) +
        facet_grid(chr ~ .id)
    dev.off()
}
foreach(i = seq(1, length(gr_list), by = 4), .packages = package_export, .export = c("txdb")) %dopar% {
    tmp <- plotAvgProf2(gr_list[i:(i + 3)],
        TxDb = txdb, upstream = 3000, downstream = 3000, ylab = "Read count frequency", xlab = "Genomic Region (5'->3')",
        conf = 0.95, resample = 1000
    ) + facet_grid(. ~ .id)
    ggsave(tmp, filename = paste0("../data/06_annotation/", names(gr_list)[i], "_avgprof", ".pdf"), width = 15, height = 10)
}
foreach(i = seq(1, length(anno_list), by = 1), .packages = package_export) %dopar% {
    pdf(paste0("../data/06_annotation/", names(anno_list)[i], "_annopie", ".pdf"))
    plotAnnoPie(anno_list[[i]])
    dev.off()
}
############################################################################################################
# annobar
############################################################################################################
get_anno_bar <- function(histone_mod_tmp, p) {
    p()
    # get elements in anno_list with names start with histone_mod_tmp
    anno_list_tmp <- anno_list[names(anno_list) %>% str_detect(histone_mod_tmp)]
    names(anno_list_tmp) <- names(anno_list_tmp) %>% str_replace(paste0(histone_mod_tmp, "_"), "")
    # order anno_list_tmp by cell_types order
    anno_list_tmp <- anno_list_tmp[cell_types]
    # plotAnnobar, empty theme, remove x axis
    sing_plot <- plotAnnoBar(anno_list_tmp, title = histone_mod_tmp) +
        scale_y_discrete(limits = cell_types) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title.x = element_blank(), panel.grid = element_blank()
        )
    sing_plot
}
with_progress({
    p <- progressor(along = histone_mod)
    annobar_list <- future_map(histone_mod, get_anno_bar, p = p) %>% setNames(histone_mod)
})
# apply rremove("ylab") to elements after the second one in annobar_list, use modify_in
annobar_list <- modify_at(
    annobar_list, seq(2, length(annobar_list)),
    function(x) x + rremove("y.text") + rremove("y.ticks") + rremove("y.title")
)
final_bar <- wrap_plots(annobar_list, ncol = 4, guides = "collect")
ggsave("a.pdf", final_bar)
foreach(
    i = seq(1, length(histone_mod), by = 1),
    .packages = package_export, .export = c("cell_types", "histone_mod", "anno_list")
) %dopar% {
    histone_mod_tmp <- histone_mod[i]
    # get elements in anno_list with names start with histone_mod_tmp
    anno_list_tmp <- anno_list[names(anno_list) %>% str_detect(histone_mod_tmp)]
    # plotAnnobar, empty theme
    annobar_list[[histone_mod_tmp]] <- plotAnnoBar(anno_list_tmp) + theme_void()
}
ggsave(tmp, filename = paste0("../data/06_annotation/", names(anno_list)[i], "_annobar", ".pdf"), width = 15, height = 10)
foreach(i = seq(1, length(anno_list), by = 1), .packages = package_export) %dopar% {
    tmp <- upsetplot(anno_list[[i]], vennpie = TRUE)
    ggsave(tmp, filename = paste0("../data/06_annotation/", names(anno_list)[i], "_annoupset", ".pdf"), width = 15, height = 10)
}
stopImplicitCluster()

library(tidyverse)
library(gridExtra)
library(ChIPseeker)
library(ggpubr)
library(ggplotify)
library(RColorBrewer)
library(patchwork)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(furrr)
library(progressr)
options(future.globals.maxSize = 891289600)
furrr_options(seed = 123)
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
cell_types <- c(
    "LT_HSC", "ST_HSC", "MPP", "CMP", "GMP", "MF",
    "GN", "Mono", "CLP", "B", "CD4", "CD8", "NK", "MEP", "EryA", "EryB"
)
# cell_types <- factor(cell_types, levels = cell_types)
# get all histone modifications in gr_list names
histone_mod <- names(gr_list) %>%
    str_extract("[:alnum:]*(?=_.*)") %>%
    unique()
# Overall annotation status
std_chr <- c(paste0("chr", 1:19), "chrM", "chrX", "chrY")
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
anno_list <- lapply(gr_list, function(x) annotatePeak(x, tssRegion = c(-3000, 3000), TxDb = txdb))
# no first exon in some gr, add a dull element to make it have the same length as others
for (i in 1:length(anno_list)) {
    # if anno_list[[i]]@annoStat$Feature does not have "1st Exon", add a dull element
    normal_order <- anno_list$H3K4me1_ST_HSC@annoStat$Feature
    if (!"1st Exon" %in% anno_list[[i]]@annoStat$Feature) {
        anno_list[[i]]@annoStat <- add_row(anno_list[[i]]@annoStat, Feature = "1st Exon", Frequency = 0.00000001)
    }
    anno_list[[i]]@annoStat <- arrange(anno_list[[i]]@annoStat, match(Feature, normal_order))
}
# Visulization
plan(multicore, workers = 20)
############################################################################################################
# cov
############################################################################################################
get_cov <- function(std_chr_tmp, p) {
    p()
    plot_line <- function(one_cell) {
        gr_list_one_cell <- gr_list[names(gr_list) %>% str_detect(paste0("_", one_cell))]
        names(gr_list_one_cell) <- names(gr_list_one_cell) %>% str_replace(paste0("_", one_cell), "")
        # add a column to all elements in gr_list_one_cell, the column is the name of the element
        line_plot <- covplot(gr_list_one_cell, ylab = one_cell, xlab = NULL, title = NULL, chrs = std_chr_tmp) +
            facet_grid(chr ~ fct_relevel(.id, sort)) + {
                if (one_cell == "LT_HSC") {
                    rremove("xy.text") + rremove("axis") + rremove("ticks") + rremove("legend") +
                        theme(strip.background.y = element_blank(), strip.text.y = element_blank(), axis.title.y = element_text(angle = 0))
                } else {
                    rremove("xy.text") + rremove("axis") + rremove("ticks") + rremove("legend") +
                        theme(
                            strip.background.y = element_blank(), strip.text.y = element_blank(),
                            strip.background.x = element_blank(), strip.text.x = element_blank(), axis.title.y = element_text(angle = 0)
                        )
                }
            }
    }
    page_plot <- map(cell_types, plot_line)
    page_plot <- wrap_plots(page_plot, ncol = 1, guides = "auto")
    page_plot
}
with_progress({
    p <- progressor(along = std_chr)
    cov_list <- future_map(std_chr, get_cov, p = p) %>% setNames(std_chr)
})
final_cov <- map(cov_list, wrap_plots, ncol = 1, guides = "collect")
# add title to each element in final_cov
final_cov <- imap(final_cov, \(x, idx) x + plot_annotation(title = idx))
pdf("../data/06_annotation_tidy/cov/cov.pdf", onefile = TRUE, height = 8)
for (i in 1:length(final_cov)) {
    print(final_cov[[i]])
}
dev.off()
############################################################################################################
# avgprof
############################################################################################################
get_avgprof_plot <- function(histone_mod_tmp, p) {
    p()
    # get elements in gr_list with names start with histone_mod_tmp
    gr_list_tmp <- gr_list[names(gr_list) %>% str_detect(histone_mod_tmp)]
    names(gr_list_tmp) <- names(gr_list_tmp) %>% str_replace(paste0(histone_mod_tmp, "_"), "")
    gr_list_tmp <- gr_list_tmp[cell_types]
    # plot avgprof
    sing_plot <- plotAvgProf2(gr_list_tmp,
        TxDb = txdb, upstream = 3000, downstream = 3000,
        conf = 0.95, resample = 1000, facet = "column",
        free_y = TRUE, xlab = NULL, ylab = histone_mod_tmp
    )
    sing_plot
}
with_progress({
    p <- progressor(along = histone_mod)
    avgprof_list <- future_map(histone_mod, get_avgprof_plot, p = p) %>% setNames(histone_mod)
})
avgprof_list <- modify(
    avgprof_list,
    function(x) x + rremove("x.text") + rremove("x.ticks") + rremove("x.title")
)
final_avgprof <- wrap_plots(avgprof_list, nrow = 4, guides = "auto")
ggsave("../data/06_annotation_tidy/avgprof/avgprof.pdf", final_avgprof, width = 14)
############################################################################################################
# annopie
############################################################################################################
# NOTE: plotAnnoPie() creates vanilla R plots, so it can't be treated in ggplot way
get_anno_pie <- function(histone_mod_tmp, p) {
    p()
    # get elements in anno_list with names start with histone_mod_tmp
    anno_list_tmp <- anno_list[names(anno_list) %>% str_detect(histone_mod_tmp %>% as.character())]
    names(anno_list_tmp) <- names(anno_list_tmp) %>% str_replace(paste0(histone_mod_tmp, "_"), "")
    # order anno_list_tmp by cell_types order
    anno_list_tmp <- anno_list_tmp[cell_types]
    # self tweak plotAnnoPie2
    plotAnnoPie2 <- function(anno, ndigit = 2, cex = 0.8, col = NA, legend.position, radius = 0.8, name, ...) {
        anno.df <- getAnnoStat(anno)
        getCols <- function(n) {
            col <- c(
                "#8dd3c7", "#ffffb3", "#bebada",
                "#fb8072", "#80b1d3", "#fdb462",
                "#b3de69", "#fccde5", "#d9d9d9",
                "#bc80bd", "#ccebc5", "#ffed6f"
            )

            col2 <- c(
                "#1f78b4", "#ffff33", "#c2a5cf",
                "#ff7f00", "#810f7c", "#a6cee3",
                "#006d2c", "#4d4d4d", "#8c510a",
                "#d73027", "#78c679", "#7f0000",
                "#41b6c4", "#e7298a", "#54278f"
            )

            col3 <- c(
                "#a6cee3", "#1f78b4", "#b2df8a",
                "#33a02c", "#fb9a99", "#e31a1c",
                "#fdbf6f", "#ff7f00", "#cab2d6",
                "#6a3d9a", "#ffff99", "#b15928"
            )
            col3[1:n]
        }
        if (is.na(col[1])) {
            col <- getCols(nrow(anno.df))
        }
        if (!all(c("Feature", "Frequency") %in% colnames(anno.df))) {
            stop("check your input...")
        }
        labels <- paste(anno.df$Feature, " (",
            round(anno.df$Frequency / sum(anno.df$Frequency) * 100, ndigit),
            "%)",
            sep = ""
        )
        p <- ggplot(anno.df, aes(x = "", y = Frequency, fill = Feature)) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar("y", start = 0) +
            theme_void() +
            theme(legend.position = "right") +
            scale_fill_manual(values = col, limits = normal_order) +
            ggtitle(name)
    }
    page_plot <- imap(anno_list_tmp, \(x, idx) plotAnnoPie2(x, name = idx))
    names(page_plot) <- names(anno_list_tmp)
    # make the order of elements in page_plot the same as cell_types
    page_plot <- page_plot[cell_types]
    page_plot <- wrap_plots(page_plot, ncol = 4, guides = "collect") + plot_annotation(title = histone_mod_tmp)
    page_plot
}
with_progress({
    p <- progressor(along = histone_mod)
    annopie_list <- map(histone_mod, get_anno_pie, p = p) %>% setNames(histone_mod)
})
# ggsave all elements in annopie_list to one pdf, each in one page of the pdf
pdf("../data/06_annotation_tidy/annopie/annopie.pdf", onefile = TRUE)
for (i in 1:length(annopie_list)) {
    print(annopie_list[[i]])
}
dev.off()
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
ggsave("../data/06_annotation_tidy/annobar/annobar.pdf", final_bar)
############################################################################################################
# heatmap
############################################################################################################
plan(multicore, workers = 20)
# install YuLab-SMU/ChIPseeker @ 15b343db0f9238dacaec12cc479438a55d996392
devtools::install_github("YuLab-SMU/ChIPseeker", ref = "15b343db0f9238dacaec12cc479438a55d996392")
# if the following code does not work, try to restart R session
unloadNamespace("ChIPseeker")
library(ChIPseeker)
get_heatmap <- function(histone_mod_tmp, p) {
    p()
    # get elements in gr_list with names start with histone_mod_tmp
    gr_list_tmp <- gr_list[names(gr_list) %>% str_detect(histone_mod_tmp)]
    names(gr_list_tmp) <- names(gr_list_tmp) %>% str_replace(paste0(histone_mod_tmp, "_"), "")
    gr_list_tmp <- gr_list_tmp[cell_types]
    # get the index of histone_mod_tmp in histone_mod
    histone_mod_idx <- which(histone_mod == histone_mod_tmp)
    col <- brewer.pal(n = 4, name = "Set1")[histone_mod_idx]
    # plot single heatmap
    if (histone_mod_idx == 1) {
        peakHeatmap(gr_list_tmp, TxDb = txdb, weightCol = "score", upstream = 3000, ylab = histone_mod_tmp, downstream = 3000, color = col)
    } else {
        peakHeatmap(gr_list_tmp, TxDb = txdb, weightCol = "score", upstream = 3000, ylab = histone_mod_tmp, title = "", downstream = 3000, color = col)
    }
}
with_progress({
    p <- progressor(along = histone_mod)
    pdf(paste0("../data/06_annotation_tidy/heatmap/", "heatmap", ".pdf"), height = 7, width = 25, onefile = TRUE)
    par(mfrow = c(4, 1))
    heatmap_list <- map(histone_mod, get_heatmap, p = p) %>% setNames(histone_mod)
    dev.off()
})
library(data.table)

source("src/column_plots.R")

# functions
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
}

# load data
tfdb.os <- readRDS("data/tfdb_os.Rds")
lmd.expressed <- readRDS("data/lmdPaper/expressedGenes/expressedGenesAll.Rds")
lmd.vst <- readRDS("data/lmdPaper/DESeq2/vst.Rds")
lmd.vstmeans <- sapply(levels(lmd.vst$stage), function(x)
  apply(
    SummarizedExperiment::assay(lmd.vst[, lmd.vst$stage == x]), 1, gm_mean))

# five accessions phenotyping
acc.pheno <- data.table(read.csv('data/Phenotype_Panicle_corrected.csv',
                               stringsAsFactors = FALSE))

# domestication genes
int.results.table <- readRDS(
  "data/fiveAccessions/deseq2/wald_tests/domestication_results_table.Rds")
dom.continent <- readRDS(
  "data/fiveAccessions/deseq2/wald_tests/domestication_by_continent_results_table.Rds")

# stage results table
stage.results.table <- readRDS(
  "data/fiveAccessions/deseq2/wald_tests/stage_species_results_table.Rds")

# 5acc expression values
tpm <- readRDS("data/fiveAccessions/tpm/tpm_with_calls.Rds")
# fix the tpm data
tpm.clean <- copy(tpm)
tpm.clean[, Species := factor(plyr::revalue(
  species,
  c(I = "O. sativa indica",
    J = "O. sativa japonica",
    R = "O. rufipogon",
    G = "O. glaberrima",
    B = "O. barthii")),
  levels = c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
             "O. barthii", "O. glaberrima"))]

acc.expressed <- readRDS("data/fiveAccessions/tpm/detected_genes.Rds")

######################
# PANICLE PHENOTYPES #
######################

# remove "Echelle" row
acc.pheno <- acc.pheno[!grep("Echelle", file_name)]

# tidy up
acc.pheno[, Accession := unlist(strsplit(
  file_name, split = "_", fixed = TRUE))[1], by = file_name]
acc.pheno[, Species := plyr::revalue(
  toupper(Accession),
  c(B88 = "O. barthii", IR64 = "O. sativa",
    NIP = "O. sativa", TOG5681 = "O. glaberrima",
    W1654 = "O. rufipogon"))]
acc.pheno[Accession == "tog5681", Accession := "Tog5681"]
acc.pheno[Accession == "Nip", Accession := "Nipponbare"]
acc.pheno[, Indiv := unlist(
  strsplit(file_name, split = "_", fixed = TRUE))[2], by = file_name]
acc.pheno[, Panicle := unlist(
  strsplit(file_name, split = "_", fixed = TRUE))[3], by = file_name]

# select data
setkey(acc.pheno, "Accession", "Indiv", "Panicle")
acc.pheno.pd.wide <- unique(acc.pheno[, .(
  Species, Accession, Indiv, Panicle,
  "Spikelets" = Sp_nb,
  "Primary branches" = SA_nb,
  "Rachis length (mm)" = PA_length * 10,
  "Secondary branches" = TA_nb)])
acc.pheno.pd <- melt(
  acc.pheno.pd.wide,
  id.vars = c("Species", "Accession", "Indiv", "Panicle"))

# order accessions and species
acc.pheno.pd[, Accession := factor(
  Accession,
  levels = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]
acc.pheno.pd[, Species := factor(
  Species,
  levels = c("O. sativa", "O. rufipogon", "O. glaberrima", "O. barthii"))]

# plot phenotypes
cols.acc <- RColorBrewer::brewer.pal(5, "Set1")[c(2, 1, 5, 4)]
acc.pheno <- ggplot(acc.pheno.pd, aes(x = Accession, y = value, fill = Species)) +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  theme_slide +
  scale_fill_manual(values = cols.acc, guide = guide_legend(title=NULL)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.text = element_text(face = "italic"),
        legend.position = "top",
        axis.text.x	= element_text(angle = 45, hjust = 1)) +
  geom_point(size = 2, alpha = 0.5, shape = 21, colour = NA,
             position = position_jitter(width = 0.4))

#######
# ap2 #
#######

# get ap2 expression data
ap2.all <- tfdb.os[Family == "AP2-EREBP", unique(Protein.ID)]
ap2 <- ap2.all[ap2.all %in% lmd.expressed]
ap2.vst <- lmd.vstmeans[ap2,]
ap2.scaled <- t(scale(t(ap2.vst), center = TRUE))

# cut into groups
hc <- hclust(dist(ap2.scaled, method = "minkowski"), method = "ward.D2")
cuts <- cutree(hc, h = 3)

# get the order for the x-axis
gene.order <- rownames(ap2.scaled)[hc$order]

# set up plot data
ap2.pd.wide <- data.table(ap2.scaled, keep.rownames = TRUE)
setnames(ap2.pd.wide, "ePBM/SBM", "ePBM/AM")
ap2.pd <- melt(ap2.pd.wide, id.vars = "rn", variable.name = "Stage",
               value.name = "Scaled reads")
ap2.pd[, symbol := oryzr::LocToGeneName(rn)$symbols, by = rn]
ap2.pd[, clust := cuts[rn]]

# set sample order
ap2.pd[, unique(Stage)]

# set symbol order
ap2.pd[, rn := factor(rn, levels = gene.order)]
ap2.pd[, order := as.numeric(rn)]
setkey(ap2.pd, order)
ap2.pd[, symbol := factor(symbol, levels = unique(symbol))]

# draw plot
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
ap2.plot <- ggplot(ap2.pd, aes(x = symbol, y = Stage, fill = `Scaled reads`)) +
  theme_slide + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_gradientn(colours = heatscale) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0),
                   limits = ap2.pd[, rev(levels(Stage))]) +
  facet_grid(. ~ clust, scales = "free_x", space = "free_x") +
  geom_raster()


#######################
# Domestication genes #
#######################

genes <- int.results.table[padj < 0.05, unique(gene)]
dom.genes <- ColumnPlot(genes) + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))

# for africa, take the top 20 in each direction
af.table <- dom.continent[padj < 0.05 & domestication == "africa"]
af.table <- af.table[order(-abs(log2FoldChange))]
af1 <- af.table[log2FoldChange > 0][1:20, unique(gene)]
af2 <- af.table[log2FoldChange < 0][1:20, unique(gene)]

# for asia, we have room to show all
as <- dom.continent[padj < 0.05 & domestication == "asia", unique(gene)]

# draw plots
as.genes <- ColumnPlot(as, species = "asian") + theme_slide +
  theme(axis.text.y = element_text(face = "italic", size = 10),
        strip.text.x = element_text(face = "italic"))
af1.genes <- ColumnPlot(af1, species = "african") + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))
af2.genes <- ColumnPlot(af2, species = "african") + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))

##################
# 5acc ap2 genes #
##################

# get ap2 expression data
ap2.acc <- ap2.all[ap2.all %in% acc.expressed]

# only genes where there is a sig de between stages in at least 1 species
ap2.acc.sig <- stage.results.table[gene %in% ap2.acc & padj < 0.5, unique(gene)]

ap2.acc.sig.plot <- ColumnPlot(ap2.acc.sig) + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))

library(ggplot2)
library(data.table)

# need the LFCs b/w PBM and SM for each species for the plot
stage.results.table <- readRDS(
  "data/fiveAccessions/deseq2/wald_tests/stage_species_results_table.Rds")
# get tpm values to map to colours
tpm.long <- readRDS("data/fiveAccessions/tpm/tpm.Rds")
ord <- c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
         "O. barthii", "O. glaberrima")
tpm.long[, accession := substr(sample, 1, 1)]
tpm.long[ , accession := factor(plyr::mapvalues(
  accession,
  from = c("R", "I", "J", "B", "G"),
  to = ord),  levels = ord)
  ]
tpm.basemean <- tpm.long[, .(mean.tpm = mean(tpm)), by = .(gene, accession)]
setkey(tpm.basemean, gene, accession)

ColumnPlot <- function(genes, species = "all"){
  
  # just for testing
  # int.results.table <- readRDS("data/fiveacc/deseq2/wald_domestication/results_table.Rds")
  # genes <- int.results.table[padj < 0.05, unique(gene)]
  
  # column plot of genes
  plot.data <- stage.results.table[gene %in% genes]
  
  # fix accession names
  ord <- c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
           "O. barthii", "O. glaberrima")
  plot.data[ , accession := factor(plyr::mapvalues(
    accession,
    from = c("rufipogon", "indica", "japonica", "barthii", "glaberrima"),
    to = ord),  levels = ord)
    ]
  
  # insert gene names if available
  plot.data[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
  plot.data[is.na(symbol), symbol := gene]
  
  # order by gene names
  # plot.data[, symbol := factor(
  #   symbol, levels = unique(symbol)[rev(order(unique(gene)))])]
  
  # dummy data to highlight wild species
  pb <- plot.data[, .(
    accession = levels(accession),
    log2FoldChange = 0,
    symbol = symbol[1])]
  pb <- pb[accession %in% c("O. rufipogon", "O. barthii")]
  pb[, accession := factor(accession, levels = accession)]
  
  pal <- RColorBrewer::brewer.pal(9, "Set1")
  heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
  
  # subset species of interest
  if (!species == "all"){
    if (species == "asian") {
      plot.data <- plot.data[accession %in% c(
        "O. rufipogon", "O. sativa indica", "O. sativa japonica")]
      pb <- pb[accession == "O. rufipogon"]
    } else if (species == "african"){
      plot.data <- plot.data[accession %in% c(
        "O. barthii", "O. glaberrima")]
      pb <- pb[accession == "O. barthii"]
    } else {
      stop("Species must be one of all, african or asian")
    }
  }
  
  # order by lfc in wild species
  soi <- plot.data[, levels(accession)[levels(accession) %in% accession][1]]
  gene.order <- plot.data[accession == soi,
                          as.character(symbol[order(log2FoldChange)])]
  plot.data[, symbol := factor(symbol, levels = rev(gene.order))]
  
  # add tpm for heatscale
  setkey(plot.data, gene, accession)
  plot.data <- tpm.basemean[plot.data]
  
  ggplot(plot.data, aes(y = symbol, x = log2FoldChange,
                        colour = log(mean.tpm + 0.5, 2))) +
    facet_grid(~ accession) +
    ylab(NULL) +
    xlab(expression(L[2]*FC["PBM–SM"] %+-% "se ("*italic(n) == "3)")) +
#    scale_colour_hue(c = 100, l = 50, h.start = 359) +
    scale_colour_gradientn(colours = heatscale,
                           limits = c(0, 10),
                           name = expression(Log[2]*TPM)) +
    geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
    geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
               size = 0.5, linetype = 2, colour = "grey") +
    geom_errorbarh(aes(xmax = log2FoldChange + lfcSE,
                       xmin = log2FoldChange - lfcSE),
                   height = 0.3, size = 0.3, colour = "black") +
    geom_point(size = 2) +
    geom_rect(data = pb, fill = NA, colour = pal[2], size = 0.5, alpha = 0.5,
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}
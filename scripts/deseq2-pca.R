##
## Modified DESeq2 plotPCA function with sample names and custom condition colors.
## Sample names will be shown next to each dot with minimal overlapping.
## The axis will display proportion of variance for each principal component.
## Tested using several DESeq2 versions from 1.12.4 to 1.14.0.
##



deseq2_pca = function(object, intgroup, ntop = 500) {

  library(DESeq2)
  library(genefilter)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)

  # PCA
  rv = rowVars(assay(object))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(object)[select, ]))

  # proportion of variance
  variance = (pca$sdev ^ 2) / (sum(pca$sdev ^ 2))
  variance = round(variance, 3) * 100

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData()")
  }

  # add the intgroup factors together to create a new grouping factor
  intgroup_df = as.data.frame(colData(object)[, intgroup, drop = FALSE])
  if (length(intgroup) > 1) {
    fac = factor(apply(intgroup_df, 1, paste, collapse = " : "))
  } else {
    fac = colData(object)[[intgroup]]
  }

  # set color sheme based on number of groups
  if (nlevels(fac) > 24) {
    colors = rainbow(nlevels(fac))
  } else if (nlevels(fac) > 2) {
    colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))
    colors = unique(colors)
    colors = head(colors, nlevels(fac))
  } else {
    colors = c("dodgerblue3", "firebrick3")
  }

  # set text font size based on number of samples (measured in mm, not points)
  num_samples = ncol(object)
  if (num_samples > 50) {
    font_size = 1.5
  } else if (num_samples > 20) {
    font_size = 2
  } else if (num_samples > 10) {
    font_size = 3
  } else {
    font_size = 4
  }

  # PCA data frame for plotting
  pca_data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = fac)

  # plot (returned, not saved)
  plot_title = paste0("PCA - ", nlevels(fac), " groups - ", num_samples, " samples")
  ggplot(data = pca_data, aes_string(x = "PC1", y = "PC2", color = "group")) +
    geom_point(size = 3) +
    scale_colour_manual(values = colors) +
    geom_text_repel(aes(label = rownames(pca_data)),
                    size = font_size, point.padding = unit(0.5, "lines"), color = "black") +
    xlab(paste0("PC1 (", variance[1], "% variance)")) +
    ylab(paste0("PC2 (", variance[2], "% variance)")) +
    ggtitle(plot_title) +
    theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())

}



# end

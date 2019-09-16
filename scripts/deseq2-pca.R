##
## Modified DESeq2 plotPCA function with sample names and custom condition colors.
## Sample names will be shown next to each dot with minimal overlapping.
## The axis will display proportion of variance for each principal component.
## Tested using several DESeq2 versions from 1.12.4 to 1.14.0.
##



deseq2_pca = function(object, intgroup, ntop = 1000) {

  suppressPackageStartupMessages({
    library(magrittr)
    library(DESeq2)
    library(genefilter)
    library(glue)
    library(RColorBrewer)
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
  })

  # run PCA
  rv = rowVars(assay(object))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(object)[select, ]))

  # determine the proportion of variance for the PCs
  pca_variance = (pca$sdev ^ 2) / (sum(pca$sdev ^ 2))
  pca_variance = round(pca_variance, 3) * 100

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
  pca_data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = fac, stringsAsFactors = FALSE)

  # randomize sample order
  pca_data = pca_data[sample(rownames(pca_data)), ]

  # plot (returned)
  ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = group), size = 4) +
    geom_text_repel(
      aes(label = rownames(pca_data)),
      size = font_size,
      point.padding = unit(0.5, "lines"),
      color = "black"
    ) +
    labs(
      title = "PCA",
      subtitle = glue("{num_samples} samples | {nlevels(fac)} groups"),
      x = glue("PC1 ({pca_variance[1]}% Variance)"),
      y = glue("PC2 ({pca_variance[2]}% Variance)")
    ) +
    scale_color_manual(values = colors) +
    theme_cowplot() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    )

}



# end

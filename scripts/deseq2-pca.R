##
## Modified DESeq2 plotPCA function with sample names and proportion of variance added.
## Sample names will be shown underneath each dot.
## The axis will display proportion of variance for each principal component.
## Tested using several DESeq2 versions from 1.2.8 to 1.12.3.
## DESeq2 plotPCA function switched from lattice to ggplot2 in version 1.5.11.
##


deseq2_pca = function(x, intgroup, ntop = 500)
{
  library(RColorBrewer)
  library(genefilter)
  library(lattice)

  # pca
  rv = rowVars(assay(x))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))

  # proportion of variance
  variance = (pca$sdev ^ 2) / (sum(pca$sdev ^ 2))
  variance = round(variance, 3) * 100

  # sample names
  names = colnames(x)

  # factor of groups
  fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))

  # colors
  if( nlevels(fac) > 24 ) {
    colors = rainbow(nlevels(fac))
  }
  else if( nlevels(fac) > 2 ) {
    colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))
    colors = unique(colors)
    colors = head(colors, nlevels(fac))
  }
  else {
    colors = c("dodgerblue3", "firebrick3")
  }

  # plot
  xyplot(
    PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), pch = 16, cex = 1.5,
    aspect = 1,
    col = colors,
    xlab = list(paste0("PC1 (", variance[1], "%)"), cex = 0.8),
    ylab = list(paste0("PC2 (", variance[2], "%)"), cex = 0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x = x, y = y, labels = names, pos = 1, offset = 0.8, cex = 0.7)
    },
    main = paste0("PCA - ", nlevels(fac), " groups - ", length(names), " samples"),
    key = list(
      space = "right",
      rect = list(col = colors),
      text = list(levels(fac)),
      rep = FALSE
    ),
    par.settings = list(par.main.text = list(cex = 0.8))
  )
}



# end

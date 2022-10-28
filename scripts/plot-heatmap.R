##
## Generate a heatmap for a matrix of values with specified genes/rows and samples/columns.
##



plot_heatmap = function(mat, row_subset, col_subset, title, col_groups = NULL, file_prefix = "heatmap") {

  suppressPackageStartupMessages({
    library(glue)
    library(pheatmap)
    library(RColorBrewer)
  })

  # confirm that the subsets exist in the full matrix
  row_subset = intersect(row_subset, rownames(mat))
  if (length(row_subset) < 2) stop("too few genes/rows")
  col_subset = intersect(col_subset, colnames(mat))
  if (length(col_subset) < 2) stop("too few samples/columns")

  # confirm that the column group labels are legitimate if provided
  if(is.null(col_groups)) {
    # pheatmap expects NA by default
    col_groups = NA
  } else {
    if (length(setdiff(col_subset, rownames(col_groups))) > 0) {
      stop("some samples/columns do not have group labels")
    }
  }

  # subset matrix to genes/rows and samples/columns of interest
  mat = mat[row_subset, col_subset]

  # heatmap cells color range: blue-white-red
  cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(50)

  # adjust row font size based on number of rows/genes
  fontsize_row = 4
  show_rownames = TRUE
  if (nrow(mat) < 60) {
    fontsize_row = 6
  }
  if (nrow(mat) > 120) {
    show_rownames = FALSE
  }

  # define group colors
  group_colors = list()
  for (g in colnames(col_groups)) {
    if (nlevels(col_groups[[g]]) > 20) {
      colors = rainbow(nlevels(col_groups[[g]]))
    } else if (nlevels(col_groups[[g]]) > 2) {
      colors = unique(c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")))
      colors = colors[1:nlevels(col_groups[[g]])]
    } else {
      colors = c("dodgerblue3", "firebrick3")
    }
    names(colors) = levels(col_groups[[g]])
    group_colors[[g]] = colors
  }

  # adjust plot width size based on number of columns/samples
  plot_width = (ncol(mat) / 4) + 3
  plot_width = round(plot_width, 1)
  plot_width = min(plot_width, 20)

  # adjust file names to reflect the data dimensions
  size_text = glue("{nrow(mat)}x{ncol(mat)}")
  filename_png = glue("{file_prefix}.{size_text}.png")
  filename_pdf = glue("{file_prefix}.{size_text}.pdf")

  # generate heatmap
  pheatmap_obj =
    pheatmap::pheatmap(
      mat, main = title, color = cell_colors, border_color = NA, scale = "row",
      fontsize_row = fontsize_row, fontsize_col = 12, show_rownames = show_rownames,
      annotation_col = col_groups, annotation_colors = group_colors, annotation_names_col = FALSE,
      cluster_rows = TRUE, cluster_cols = TRUE,
      treeheight_row = 25, treeheight_col = 25, silent = TRUE
    )

  # save heatmap as png
  message("generate heatmap: ", filename_png)
  png(filename_png, width = plot_width, height = 8, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(pheatmap_obj$gtable)
  dev.off()
  Sys.sleep(1)

  # save heatmap as pdf
  message("generate heatmap: ", filename_pdf)
  pdf(filename_pdf, width = plot_width, height = 8, version = "1.7")
  grid::grid.newpage()
  grid::grid.draw(pheatmap_obj$gtable)
  dev.off()
  Sys.sleep(1)

}



# end

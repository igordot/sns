##
## Generate a heatmap for a matrix of values with specified genes/rows and samples/columns.
##



plot_heatmap = function(mat, row_subset, col_subset, title, file_prefix = "plot.heatmap") {

  library(glue)
  library(pheatmap)

  # confirm that the subsets exist in the full matrix
  row_subset = intersect(row_subset, rownames(mat))
  if (length(row_subset) < 2) stop("too few genes/rows")
  col_subset = intersect(col_subset, colnames(mat))
  if (length(col_subset) < 2) stop("too few samples/columns")

  # subset matrix to genes/rows and samples/columns of interest
  mat = mat[row_subset, col_subset]

  # labels
  size_text = glue("{nrow(mat)}x{ncol(mat)}")
  title = glue("{title} - {size_text}")
  filename_png = glue("{file_prefix}.{size_text}.png")
  filename_pdf = glue("{file_prefix}.{size_text}.pdf")

  # heatmap cells color range: blue-white-red
  cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(50)

  # adjust row font size based on number of rows
  fontsize_row = 4
  show_rownames = TRUE
  if (nrow(mat) < 80) {
    fontsize_row = 6
  }
  if (nrow(mat) > 150) {
    show_rownames = FALSE
  }

  # png heatmap with clustering
  message("generate heatmap: ", filename_png)
  pheatmap(mat, color = cell_colors, border_color = NA, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
           main = title, fontsize_row = fontsize_row, fontsize_col = 12, show_rownames = show_rownames,
           filename = filename_png, width = 8, height = 8)
  Sys.sleep(1)

  # pdf heatmap with clustering
  message("generate heatmap: ", filename_pdf)
  pheatmap(mat, color = cell_colors, border_color = NA, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
           main = title, fontsize_row = fontsize_row, fontsize_col = 12, show_rownames = show_rownames,
           filename = filename_pdf, width = 8, height = 8)
  Sys.sleep(1)

}



# end

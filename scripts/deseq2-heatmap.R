##
## Generate a heatmap for a matrix of values with specified genes/rows and samples/columns.
##



deseq2_heatmap = function(mat, genes, samples, title, file_suffix) {

  library(pheatmap)

  # subset matrix to genes and samples of interest
  mat = mat[genes, samples]

  # labels
  size_text = paste0(length(genes), "x", length(samples))
  title = paste0(title, " - ", size_text)
  filename_png = paste0("plot.heatmap.", file_suffix, ".", size_text, ".png")
  filename_pdf = paste0("plot.heatmap.", file_suffix, ".", size_text, ".pdf")

  # heatmap cells color range: blue-white-red
  cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(50)

  # adjust row font size based on gene list size
  fontsize_row = 4
  show_rownames = TRUE
  if (length(genes) < 80) {
    fontsize_row = 6
  }
  if (length(genes) > 150) {
    show_rownames = FALSE
  }

  # heatmap with clustering
  message("generate heatmap: ", filename_png)
  pheatmap(mat, color = cell_colors, border_color = NA, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",
  main = title, fontsize_row = fontsize_row, fontsize_col = 12, show_rownames = show_rownames,
  filename = filename_png, width = 12, height = 8)
  Sys.sleep(1)

  # heatmap with clustering
  message("generate heatmap: ", filename_pdf)
  pheatmap(mat, color = cell_colors, border_color = NA, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",
  main = title, fontsize_row = fontsize_row, fontsize_col = 12, show_rownames = show_rownames,
  filename = filename_pdf, width = 12, height = 8)
  Sys.sleep(1)

}



# end

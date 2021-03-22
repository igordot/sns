##
## Generate a volcano plot from a data frame of genes, fold changes (log-scale), and p-values.
##



plot_volcano = function(stats_df,
                        gene_col,
                        fc_col,
                        p_col,
                        fc_cutoff = 0,
                        p_cutoff = 0.05,
                        fc_label = "Fold Change",
                        p_label = "P-Value",
                        title = "Volcano Plot",
                        n_top_genes = 10,
                        file_prefix = "volcano") {

  suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(glue)
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
  })

  if (nrow(stats_df) < 100) stop("stats table too short")

  # set up the data frame for the volcano plot and create a column to define the significant genes
  volcano_tbl =
    stats_df %>%
    dplyr::select(gene = .data[[gene_col]], fc = .data[[fc_col]], p = .data[[p_col]]) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(sig = dplyr::if_else(abs(fc) > fc_cutoff & p < p_cutoff, "yes", "no"))

  # plot
  volcano_plot =
    ggplot(volcano_tbl, aes(x = fc, y = -log10(p), color = sig)) +
    geom_point(size = 1.2) +
    geom_text_repel(
      data = head(volcano_tbl, n_top_genes),
      aes(label = gene),
      color = "black",
      size = 4,
      point.padding = 0.1
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = title,
      x = fc_label,
      y = p_label
    ) +
    scale_color_manual(values = c("slategrey", "firebrick2")) +
    theme_cowplot() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )

  volcano_png = glue("{file_prefix}.png")
  message("generate volcano plot: ", volcano_png)
  save_plot(filename = volcano_png, plot = volcano_plot, base_height = 8, base_width = 6)
  Sys.sleep(1)

  volcano_pdf = glue("{file_prefix}.pdf")
  message("generate volcano plot: ", volcano_pdf)
  save_plot(filename = volcano_pdf, plot = volcano_plot, base_height = 8, base_width = 6)
  Sys.sleep(1)

}



# end

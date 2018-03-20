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
                        p_label = "P Value",
                        title = "Volcano Plot",
                        n_top_genes = 10,
                        file_prefix = "plot.volcano") {

  library(glue)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)

  if (nrow(stats_df) < 100) stop("stats table too short")

  # set up the data frame for the volcano plot
  volcano_df = data.frame(gene = stats_df[, gene_col],
                          fc = stats_df[, fc_col],
                          p = stats_df[, p_col],
                          stringsAsFactors = FALSE)

  # create a column to define the significant genes (by default, p<0.05 and no fold change cutoff)
  volcano_df$sig = as.factor(abs(volcano_df$fc) > fc_cutoff & volcano_df$p < p_cutoff)

  # plot
  volcano_plot = ggplot(volcano_df, aes(x = fc, y = -log10(p), color = sig)) +
    geom_point(size = 1.2) +
    xlab(fc_label) +
    ylab(p_label) +
    theme(legend.position = "none") +
    geom_text_repel(data = head(volcano_df, n_top_genes), aes(label = gene),
                    color = "black", size = 4, point.padding = 0.1) +
    ggtitle(title) +
    scale_color_manual(values = c("slategrey", "firebrick2"))

  volcano_png = glue("{file_prefix}.png")
  message("generate volcano plot: ", volcano_png)
  save_plot(volcano_png, volcano_plot, base_height = 8, base_width = 6)
  Sys.sleep(1)

  volcano_pdf = glue("{file_prefix}.pdf")
  message("generate volcano plot: ", volcano_pdf)
  save_plot(volcano_pdf, volcano_plot, base_height = 8, base_width = 6)
  Sys.sleep(1)

}



# end

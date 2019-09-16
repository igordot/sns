##
## Calculate results and export them in various formats with proper names.
##



deseq2_compare = function(deseq_dataset, contrast = NULL, name = NULL) {

  suppressPackageStartupMessages({
    library(magrittr)
    library(glue)
    library(DESeq2)
    library(ashr)
    library(readr)
    library(writexl)
  })

  # create separate directories for certain output files
  r_dir = "r-data"
  if (!dir.exists(r_dir)) dir.create(r_dir)
  heatmaps_dir = "heatmaps"
  if (!dir.exists(heatmaps_dir)) dir.create(heatmaps_dir)
  volcano_dir = "volcano-plots"
  if (!dir.exists(volcano_dir)) dir.create(volcano_dir)

  # calculate results (using contrast or name, depending on what is given)
  # since v1.16 (11/2016), lfcShrink function performs fold change shrinkage and addMLE is for backward compatibility
  # since v1.18 (11/2017), two alternative shrinkage estimators added ("apeglm" and "ashr")
  # "normal" and "ashr" can be used with arbitrary specified contrast ("apeglm" does not shrink multiple coefficients)
  # for "ashr", if res is provided, then coef and contrast are ignored
  if(!is.null(contrast)) {
    res_unshrunk = results(deseq_dataset, contrast = contrast, cooksCutoff = FALSE)
    res = lfcShrink(deseq_dataset, contrast = contrast, res = res_unshrunk, type = "ashr")
    # extract results name
    pattern = paste(".*", contrast[1], " ", sep = "")
    res_name = gsub(pattern = pattern, replacement = "", x = mcols(res)[2,2])
    samples_comp = rownames(subset(deseq_dataset@colData, group %in% contrast[2:3]))
  } else {
    # not tested in combination with lfcShrink
    res = results(deseq_dataset, name = name, cooksCutoff = FALSE, addMLE = TRUE)
    res_name = name
  }

  # file suffix based on comparison name
  file_suffix = gsub(pattern = " ", replacement = "-", x = res_name)

  # sort results so most significant are first
  res = res[order(res$padj, res$pvalue, -res$baseMean), ]

  # save unmodified results object
  res_rds = glue("{r_dir}/deseq2.res.{file_suffix}.rds")
  saveRDS(res, file = res_rds)
  message("save results object: ", res_rds)
  Sys.sleep(1)

  # save the unmodified results table as csv
  res_tbl = as_tibble(res, rownames = "gene") %>% dplyr::arrange(padj, pvalue, -baseMean)
  res_csv = glue("dge.{file_suffix}.csv")
  write_excel_csv(res_tbl, path = res_csv)
  message("save results csv: ", res_csv)
  Sys.sleep(1)

  # add unshrunk fold change to results
  res_unshrunk_tbl = as_tibble(res_unshrunk, rownames = "gene")
  res_unshrunk_tbl = dplyr::select(res_unshrunk_tbl, gene, log2FCunshrunk = log2FoldChange)
  res_tbl = left_join(res_tbl, res_unshrunk_tbl, by = "gene")

  # format results for excel export
  res_tbl =
    res_tbl %>%
    dplyr::mutate(
      baseMean       = round(baseMean, 1),
      log2FC         = round(log2FoldChange, 2),
      log2FCunshrunk = round(log2FCunshrunk, 2),
      pvalue         = if_else(pvalue < 0.00001, pvalue, round(pvalue, 5)),
      padj           = if_else(padj < 0.00001, padj, round(padj, 5))
    ) %>%
    dplyr::select(gene, baseMean, log2FC, log2FCunshrunk, pvalue, padj)

  message("num genes padj<0.90: ", nrow(subset(res_tbl, padj < 0.9)))
  message("num genes padj<0.20: ", nrow(subset(res_tbl, padj < 0.2)))
  message("num genes padj<0.05: ", nrow(subset(res_tbl, padj < 0.05)))
  message("num genes padj<0.01: ", nrow(subset(res_tbl, padj < 0.01)))

  # save differential expression results in Excel format
  res_xlsx = paste0("dge.", gsub(pattern = " ", replacement = "-", x = res_name), ".xlsx")
  write_xlsx(setNames(list(res_tbl), strtrim(res_name, 31)), path = res_xlsx)
  message("results genes: ", nrow(res_tbl))
  message("save results xlsx: ", res_xlsx)
  Sys.sleep(1)

  # save significant (padj<0.05) differential expression results in Excel format
  res_padj005_xlsx = gsub(pattern = ".xlsx", replacement = ".q005.xlsx", x = res_xlsx)
  res_padj005_df = subset(res_tbl, padj < 0.05)
  write_xlsx(setNames(list(res_padj005_df), strtrim(res_name, 31)), path = res_padj005_xlsx)
  message("save filtered results xlsx: ", res_padj005_xlsx)
  Sys.sleep(1)

  # generate volcano plot
  plot_volcano(stats_df = res_tbl, gene_col = "gene", fc_col = "log2FC", p_col = "padj", p_cutoff = 0.05,
               fc_label = "Fold Change (log2)", p_label = "Adjusted P Value (-log10)",
               title = res_name,
               file_prefix = glue("{volcano_dir}/plot.volcano.{file_suffix}"))

  # heatmap variance stabilized values matrix
  vsd = assay(varianceStabilizingTransformation(deseq_dataset, blind = TRUE))

  # all samples and the subset used for the comparison
  samples_all = colnames(deseq_dataset)
  samples_comp = samples_all
  if(!is.null(contrast)) samples_comp = rownames(subset(deseq_dataset@colData, group %in% contrast[2:3]))

  # heatmap gene subsets (list with genes, plot title, and file suffix)
  hmg = list()
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% head(50) %>% pull(gene),
                                title = "50 Most Significant",
                                file_suffix = "top")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% head(100) %>% pull(gene),
                                title = "100 Most Significant",
                                file_suffix = "top")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% head(1000) %>% pull(gene),
                                title = "1000 Most Significant",
                                file_suffix = "top")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% dplyr::filter(padj < 0.10) %>% pull(gene),
                                title = "q < 0.1",
                                file_suffix = "q010")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% dplyr::filter(padj < 0.05) %>% pull(gene),
                                title = "q < 0.05",
                                file_suffix = "q005")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% dplyr::filter(padj < 0.01) %>% pull(gene),
                                title = "q < 0.01",
                                file_suffix = "q001")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% dplyr::filter(padj < 0.001) %>% pull(gene),
                                title = "q < 0.001",
                                file_suffix = "q0001")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% dplyr::filter(pvalue < 0.05) %>% pull(gene),
                                title = "p < 0.05",
                                file_suffix = "p005")
  hmg[[length(hmg) + 1]] = list(genes = res_tbl %>% dplyr::filter(pvalue < 0.01) %>% pull(gene),
                                title = "p < 0.01",
                                file_suffix = "p001")

  # process every gene subset
  for (i in 1:length(hmg)) {

    # generate title and file suffix
    hm_title = glue("{res_name} - {hmg[[i]]$title}")
    hm_file_prefix = glue("{heatmaps_dir}/plot.heatmap.{file_suffix}.{hmg[[i]]$file_suffix}")

    # generate heatmaps if gene list is not too small or big
    if (length(hmg[[i]]$genes) > 10 && length(hmg[[i]]$genes) < 3000) {

      # generate heatmap using all samples
      plot_heatmap(mat = vsd, row_subset = hmg[[i]]$genes, col_subset = samples_all,
        title = hm_title, file_prefix = hm_file_prefix)

      # generate heatmap using a subset of samples used for the comparison
      if (length(samples_comp) < length(samples_all)) {
        plot_heatmap(mat = vsd, row_subset = hmg[[i]]$genes, col_subset = samples_comp,
          title = hm_title, file_prefix = hm_file_prefix)
      }

    }

  }

}



# end

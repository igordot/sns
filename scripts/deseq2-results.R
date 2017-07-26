##
## Calculate results and export them in various formats with proper names.
##



deseq2_results = function(deseq_dataset, contrast = NULL, name = NULL) {

  library(DESeq2)
  library(xlsx)

  # calculate results (using contrast or name, depending on what is given)
  if(!is.null(contrast)) {
    res = results(deseq_dataset, contrast = contrast, cooksCutoff = FALSE, addMLE = TRUE)
    # extract results name
    pattern = paste(".*", contrast[1], " ", sep = "")
    res_name = gsub(pattern = pattern, replacement = "", x = mcols(res)[2,2])
    samples_comp = rownames(subset(deseq_dataset@colData, group %in% contrast[2:3]))
  } else {
    res = results(deseq_dataset, name = name, cooksCutoff = FALSE, addMLE = TRUE)
    res_name = name
  }

  # file suffix based on comparison name
  file_suffix = gsub(pattern = " ", replacement = "-", x = res_name)

  # sort results so most significant are first
  res = res[order(res$padj, res$pvalue, -res$baseMean), ]

  # save results object
  res_rds = paste0("deseq2.res.", file_suffix, ".rds")
  saveRDS(dds, file = res_rds)
  message("save results object: ", res_rds)

  # save results as csv
  res_csv = paste0("dge.", file_suffix, ".csv")
  write.csv(as.data.frame(res), file = res_csv)
  message("save results csv: ", res_csv)

  # format results for excel export
  res_df = as.data.frame(res)
  res_df$gene           = rownames(res_df)
  res_df$baseMean       = round(res_df$baseMean, 1)
  res_df$log2FoldChange = round(res_df$log2FoldChange, 5)
  res_df$lfcMLE         = round(res_df$lfcMLE, 5)
  res_df$pvalue         = round(res_df$pvalue, 15)
  res_df$padj           = round(res_df$padj, 15)
  res_df = subset(res_df, select = c("gene", "baseMean", "log2FoldChange", "lfcMLE", "pvalue", "padj"))

  # rename columns
  names(res_df)[names(res_df) == "log2FoldChange"] = "log2FC"
  names(res_df)[names(res_df) == "lfcMLE"] = "log2FCunshrunk"

  # save results as xlsx
  res_xlsx = paste0("dge.", gsub(pattern = " ", replacement = "-", x = res_name), ".xlsx")
  write.xlsx2(x = res_df, file = res_xlsx, sheetName = res_name, col.names = TRUE, row.names = FALSE)
  message("results genes: ", nrow(res_df))
  message("save results xlsx: ", res_xlsx)

  # save results as xlsx (padj<0.05)
  res_padj005_xlsx = gsub(pattern = ".xlsx", replacement = ".q005.xlsx", x = res_xlsx)
  res_padj005_df = subset(res_df, padj < 0.05)
  write.xlsx2(x = res_padj005_df, file = res_padj005_xlsx, sheetName = res_name, col.names = TRUE, row.names = FALSE)
  message("num genes padj<0.9:  ", nrow(subset(res_df, padj < 0.9)))
  message("num genes padj<0.2:  ", nrow(subset(res_df, padj < 0.2)))
  message("num genes padj<0.05: ", nrow(subset(res_df, padj < 0.05)))
  message("num genes padj<0.01: ", nrow(subset(res_df, padj < 0.01)))
  message("save filtered results xlsx: ", res_padj005_xlsx)

  # heatmap values matrix
  vsd = assay(varianceStabilizingTransformation(dds, blind = TRUE))

  # all samples and the subset used for the comparison
  samples_all = colnames(deseq_dataset)
  samples_comp = samples_all
  if(!is.null(contrast)) samples_comp = rownames(subset(deseq_dataset@colData, group %in% contrast[2:3]))

  # heatmap gene subsets (list with genes, plot title, and file suffix)
  hmg = list()
  hmg[[length(hmg) + 1]] = list(genes = rownames(res_df)[1:50], title = "50 Most Significant", file_suffix = "top")
  hmg[[length(hmg) + 1]] = list(genes = rownames(res_df)[1:100], title = "100 Most Significant", file_suffix = "top")
  hmg[[length(hmg) + 1]] = list(genes = rownames(res_df)[1:1000], title = "1000 Most Significant", file_suffix = "top")
  hmg[[length(hmg) + 1]] = list(genes = rownames(subset(res_df, padj < 0.10)), title = "q < 0.1", file_suffix = "q010")
  hmg[[length(hmg) + 1]] = list(genes = rownames(subset(res_df, padj < 0.05)), title = "q < 0.05", file_suffix = "q005")
  hmg[[length(hmg) + 1]] = list(genes = rownames(subset(res_df, padj < 0.01)), title = "q < 0.01", file_suffix = "q001")
  hmg[[length(hmg) + 1]] = list(genes = rownames(subset(res_df, padj < 0.001)), title = "q < 0.001", file_suffix = "q0001")
  hmg[[length(hmg) + 1]] = list(genes = rownames(subset(res_df, pvalue < 0.05)), title = "p < 0.05", file_suffix = "p005")
  hmg[[length(hmg) + 1]] = list(genes = rownames(subset(res_df, pvalue < 0.01)), title = "p < 0.01", file_suffix = "p001")

  # process every gene subset
  for (i in 1:length(hmg)) {

    # generate title and file suffix
    hm_title = paste0(res_name, " - ", hmg[[i]]$title)
    hm_file_suffix = paste0(file_suffix, ".", hmg[[i]]$file_suffix)

    # generate heatmaps if gene list is not too small or big
    if (length(hmg[[i]]$genes) > 5 && length(hmg[[i]]$genes) < 3000) {

      # generate heatmap using all samples
      deseq2_heatmap(mat = vsd, genes = hmg[[i]]$genes, samples = samples_all, title = hm_title, file_suffix = hm_file_suffix)

      # generate heatmap using a subset of samples used for the comparison
      if (length(samples_comp) < length(samples_all)) {
        deseq2_heatmap(mat = vsd, genes = hmg[[i]]$genes, samples = samples_comp, title = hm_title, file_suffix = hm_file_suffix)
      }

    }

  }

}



# end

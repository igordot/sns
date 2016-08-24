##
## Calculate results and export them in various formats with proper names.
##


deseq2_results = function(DESeqDataSet, contrast=NULL, name=NULL)
{
  library(xlsx)

  # calculate results (using contrast or name, depending on what is given)
  if(!is.null(contrast)) {
    res = results(DESeqDataSet, contrast=contrast, cooksCutoff=FALSE, addMLE=TRUE)
    # extract results name
    pattern = paste(".*", contrast[1], " ", sep="")
    res.name = gsub(pattern=pattern, replacement="", x=mcols(res)[2,2])
  }
  else {
    res = results(DESeqDataSet, name=name, cooksCutoff=FALSE, addMLE=TRUE)
    res.name = name
  }

  # sort results
  res = res[order(res$padj, res$pvalue, -res$baseMean),]

  # save results object
  res.rdata = paste("res.", gsub(pattern=" ", replacement="-", x=res.name), ".RData", sep="")
  save(res, file=res.rdata)
  message("save results object: ", res.rdata)

  # save results as csv
  res.csv = paste0("dge.", gsub(pattern=" ", replacement="-", x=res.name), ".csv")
  write.csv(as.data.frame(res), file=res.csv)
  message("save results csv: ", res.csv)

  # format results for excel export
  res.df = as.data.frame(res)
  res.df$gene           = rownames(res.df)
  res.df$baseMean       = round(res.df$baseMean, 1)
  res.df$log2FoldChange = round(res.df$log2FoldChange, 5)
  res.df$lfcMLE         = round(res.df$lfcMLE, 5)
  res.df$pvalue         = round(res.df$pvalue, 15)
  res.df$padj           = round(res.df$padj, 15)
  res.df = subset(res.df, select = c("gene", "baseMean", "log2FoldChange", "lfcMLE", "pvalue", "padj"))

  # rename columns
  names(res.df)[names(res.df)=="log2FoldChange"] = "log2FC"
  names(res.df)[names(res.df)=="lfcMLE"] = "log2FCunshrunk"

  # save results as xlsx
  res.xlsx = paste0("dge.", gsub(pattern=" ", replacement="-", x=res.name), ".xlsx")
  write.xlsx2(x=res.df, file=res.xlsx, sheetName=res.name, col.names=TRUE, row.names=FALSE)
  message("results genes: ", nrow(res.df))
  message("save results xlsx: ", res.xlsx)

  # save results as xlsx (padj<0.05)
  res.padj005.xlsx = gsub(pattern=".xlsx", replacement=".padj005.xlsx", x=res.xlsx)
  write.xlsx2(x=subset(res.df, padj<0.05), file=res.padj005.xlsx, sheetName=res.name, col.names=TRUE, row.names=FALSE)
  message("num genes padj<0.9:  ", nrow(subset(res.df, padj<0.9)))
  message("num genes padj<0.2:  ", nrow(subset(res.df, padj<0.2)))
  message("num genes padj<0.05: ", nrow(subset(res.df, padj<0.05)))
  message("num genes padj<0.01: ", nrow(subset(res.df, padj<0.01)))
  message("save filtered results xlsx: ", res.padj005.xlsx)
}



# end

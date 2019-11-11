##
## Perform gene set enrichment analysis from a data frame of genes and ranks.
##



gse_fgsea = function(stats_df, gene_col, rank_col, species, title = "", file_prefix = "gse") {

  suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(stringr)
    library(glue)
    library(readr)
    library(msigdbr)
    library(fgsea)
    library(ggplot2)
    library(cowplot)
  })

  if (nrow(stats_df) < 100) stop("ranks table too short")

  # if genome build was specified, convert to msigdbr-compatible species name
  if (species %in% c("hg19", "hg38")) {
    species = "Homo sapiens"
  } else if (species %in% c("mm9", "mm10")) {
    species = "Mus musculus"
  } else {
    # species = NULL
  }

  # extract species-specific gene sets
  genesets_tbl = msigdbr(species = species) %>% dplyr::mutate(gs_name = str_trunc(gs_name, 100))

  # check for sufficient number of gene set genes (if an unsupported species was used)
  genes_msigdb = genesets_tbl %>% dplyr::pull(gene_symbol) %>% unique()
  if (length(genes_msigdb) < 1000) {
    warning("gene sets table too short")
    return(NULL)
  }

  # set up the gene ranks data frame
  ranks_tbl =
    stats_df %>%
    dplyr::select(gene = gene_col, rank = rank_col) %>%
    tidyr::drop_na() %>%
    dplyr::filter(abs(rank) > 0) %>%
    dplyr::distinct() %>%
    dplyr::arrange(-rank)

  if (length(unique(ranks_tbl$gene)) < length(ranks_tbl$gene)) stop("ranks table contains duplicate genes")
  if (nrow(ranks_tbl) < 100) stop("ranks table too short after removing 0s")

  # save the ranks table
  write_excel_csv(ranks_tbl, path = glue("{file_prefix}.ranks.csv"))
  Sys.sleep(1)

  # convert the gene ranks data frame to a vector
  ranks = tibble::deframe(ranks_tbl)

  # check that there is a substantial overlap between genes in the gene set and ranks tables
  if (length(intersect(genes_msigdb, names(ranks))) < min(length(genes_msigdb), length(ranks)) / 2) {
    message("num ranked genes: ", length(ranks))
    message("num gene set genes: ", length(genes_msigdb))
    message("num overlap genes: ", length(intersect(genes_msigdb, names(ranks))))
    stop("most ranked genes are not represented in the gene set genes")
  }

  # specify categories of interest (split C2 and C5 by sub-categories)
  geneset_cats =
    c(
      "Hallmark" = "H",
      "Chemical and Genetic Perturbations" = "CGP",
      "KEGG" = "CP:KEGG",
      "Pathway Interaction Database" = "CP:PID",
      "Reactome" = "CP:REACTOME",
      "GO Biological Process" = "BP",
      "GO Cellular Component" = "CC",
      "GO Molecular Function" = "MF",
      "Oncogenic" = "C6"
    )

  # run gene set enrichment for each category
  for (geneset_name in names(geneset_cats)) {

    message("geneset enrichment: ", geneset_name)

    # define output file names
    geneset_cat = geneset_cats[geneset_name]
    geneset_cat_str = str_replace(geneset_cat, ":", "-")
    geneset_prefix = glue("{file_prefix}.fgsea.{geneset_cat_str}")

    # extract relevant gene sets
    geneset_list = genesets_tbl %>% filter(gs_cat == geneset_cat)
    if (nrow(geneset_list) == 0) {
      geneset_list = genesets_tbl %>% filter(gs_subcat == geneset_cat)
    }
    geneset_list = split(x = geneset_list$gene_symbol, f = geneset_list$gs_name)

    # run preranked gene set enrichment analysis (using tryCatch to ignore potential errors)
    bpparam = BiocParallel::MulticoreParam(4)
    set.seed(99)
    tryCatch(
      {
        # fgsea_res = fgsea(pathways = geneset_list, stats = ranks, minSize = 10, maxSize = 500, nperm = 1e5, BPPARAM = bpparam)
        fgsea_res = fgseaMultilevel(pathways = geneset_list, stats = ranks, minSize = 10, BPPARAM = bpparam)
      },
      error = function(e) {
        message("fgseaMultilevel error:", conditionMessage(e))
      }
    )

    # export fgsea results table (fgseaMultilevel replaces nMoreExtreme with log2err)
    fgsea_res$leading_edge_genes = sapply(fgsea_res$leadingEdge, paste, collapse = "|")
    fgsea_tbl =
      fgsea_res %>%
      dplyr::filter(pval < 0.1) %>%
      dplyr::arrange(padj, pval, -abs(NES)) %>%
      dplyr::mutate(
        pval = if_else(pval < 0.00001, pval, round(pval, 5)),
        padj = if_else(padj < 0.00001, padj, round(padj, 5)),
        ES = round(ES, 5),
        NES = round(NES, 5)
      ) %>%
      dplyr::select(-log2err, -leadingEdge)
    write_excel_csv(fgsea_tbl, path = glue("{geneset_prefix}.csv"))
    Sys.sleep(1)

    # filter fgsea results table for plotting
    fgsea_top_tbl =
      fgsea_tbl %>%
      tidyr::drop_na() %>%
      dplyr::filter(padj < 0.2) %>%
      dplyr::arrange(padj, pval, -abs(NES)) %>%
      dplyr::mutate(nes_dir = if_else(NES > 0, "Pos", "Neg")) %>%
      head(50) %>%
      dplyr::arrange(-NES) %>%
      dplyr::mutate(pathway = str_trunc(pathway, 50))

    # generate bar plot if any pathways remain after filtering
    if (nrow(fgsea_top_tbl) > 1) {
      fgsea_plot =
        ggplot(fgsea_top_tbl, aes(x = reorder(pathway, NES), y = NES)) +
        geom_col(aes(fill = nes_dir, alpha = padj)) +
        geom_hline(yintercept = 0) +
        labs(
          title = glue("{title}"),
          subtitle = glue("MSigDB {geneset_name} Gene Sets\n50 Most Significant"),
          x = "Gene Set",
          y = "Normalized Enrichment Score"
        )  +
        theme_cowplot() +
        theme(
          aspect.ratio = 3,
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        ) +
        coord_flip() +
        scale_fill_manual(name = "NES", values = c(Neg = "#053061", Pos = "#E41A1C")) +
        scale_alpha_continuous(name = "Adj P-Val", range = c(1, 0.2))
      save_plot(filename = glue("{geneset_prefix}.barplot.png"), plot = fgsea_plot, base_height = 10, base_width = 12)
      Sys.sleep(1)
      save_plot(filename = glue("{geneset_prefix}.barplot.pdf"), plot = fgsea_plot, base_height = 10, base_width = 12)
      Sys.sleep(1)
    }

  }

}



# end

library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

path_to_counts ="./data/Counts/subread_counts.txt"

# NK list 
#  https://aacrjournals.org/cancerimmunolres/article/7/7/1162/469488/A-Gene-Signature-Predicting-Natural-Killer-Cell
NK_signature = c("CD160","CD244","CTSW","FASLG","GZMA","GZMB",
                 "GZMH","IL18RAP","IL2RB",'KIR2DL4',"KLRB1",
                 "KLRC3",'KLRD1',"KLRF1","KLRK1","NCR1","NKG7",
                 "PRF1","XCL1","XCL2")

# genes which we want to have low expression
not_NK = c("CD19","CD8B","CD4","CD14","MS4A1") # MS4A1 = cd20

ggplot2::theme_set(theme_bw()+
                     theme(axis.line.x.top = element_blank(),
                           axis.line.y.right = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           axis.line.x.bottom = element_line(),
                           axis.line.y.left = element_line(),
                           legend.key.height = unit(1, "mm")))

get_DEGs = function(dgelrt, logFC.thresh = 2, p.thresh = 0.05){
  tmp = as.data.frame(dgelrt) %>%
    mutate(FDR = p.adjust(PValue, "BH")) %>%
    mutate(result = ifelse(logFC < -logFC.thresh & FDR < p.thresh,
                           "DOWN",
                           ifelse(logFC > logFC.thresh & FDR < p.thresh,
                                  "UP", "NS"))) %>%
    arrange(desc(abs(logFC)))
  up = tmp$genes[tmp$result == "UP"]
  down = tmp$genes[tmp$result == "DOWN"]
  return(list("up" = up, "down" = down))
}

get_num_DEGs = function(dgelrt, logFC.thresh = 2, p.thresh = 0.05){
  as.data.frame(dgelrt) %>%
    mutate(FDR = p.adjust(PValue, "BH")) %>%
    mutate(result = ifelse(logFC < -logFC.thresh & FDR < p.thresh,
                           "DOWN",
                           ifelse(logFC > logFC.thresh & FDR < p.thresh,
                                  "UP", "NS"))) %>%
    pull(result) %>%
    table
}

plot_MD = function(dgelrt, title = "", x_pct = 0.98, y_pct = 0.99, up_offset = 0.5,
                   abs_lfc_thresh = 0){
  # plots an MD plot using the results from an edgeR QLFtest.
  # x_pct and y_pct specify at which quantiles the # DE genes should be plotted
  # title is the title of the plot

  tbl = dgelrt$table %>%
    dplyr::mutate(FDR = p.adjust(PValue, method = "BH"))  %>%
    dplyr::mutate(gene = dgelrt$genes$genes)
  n_downreg = sum(tbl$FDR < 0.05 & tbl$logFC < 0 & abs(tbl$logFC) > abs_lfc_thresh)
  n_upreg = sum(tbl$FDR < 0.05 & tbl$logFC > 0 & abs(tbl$logFC) > abs_lfc_thresh)
  x_loc = quantile(tbl$logCPM, x_pct) %>% unname
  y_loc = quantile(tbl$logFC, y_pct) %>% unname
  tbl %>%
    mutate(`contrast significant` = FDR < 0.05 & abs(logFC) > abs_lfc_thresh) %>%
    ggplot(aes(x = logCPM, y = logFC, 
               size = `contrast significant`,
               color = `contrast significant`,
               label = gene))+
    annotate("text", x = x_loc, y = y_loc, hjust = 0,
             label = paste("# down =", n_downreg))+
    annotate("text", x = x_loc, y = y_loc+up_offset, hjust = 0,
             label = paste("# up =", n_upreg))+
    scale_size_manual(values = c(0.1,1))+
    scale_color_manual(values = c("black", "blue"))+
    geom_point()+
    labs(title = title)
}
plot_volcano = function(dgelrt, title = "", x_pct = 0.98, y_pct = 0.99, up_offset = 0.5,
                        abs_lfc_thresh = 0){
  # plots a volcano plot using the results from an edgeR QLFtest.
  # x_pct and y_pct specify at which quantiles the # DE genes should be plotted
  # title is the title of the plot
  
  tbl = dgelrt$table %>%
    mutate(FDR = p.adjust(PValue, method = "BH"))  %>%
    mutate(neglog10FDR = -log10(FDR))
  n_downreg = sum(tbl$FDR < 0.05 & tbl$logFC < 0 & abs(tbl$logFC) > abs_lfc_thresh)
  n_upreg = sum(tbl$FDR < 0.05 & tbl$logFC > 0 & abs(tbl$logFC) > abs_lfc_thresh)
  x_loc = quantile(tbl$logFC, x_pct) %>% unname
  y_loc = quantile(tbl$neglog10FDR, y_pct) %>% unname
  tbl %>%
    mutate(`contrast significant` = FDR < 0.05 & abs(logFC) > abs_lfc_thresh) %>%
    ggplot(aes(x = logFC, y = neglog10FDR, 
               size = `contrast significant`,
               color = `contrast significant`))+
    annotate("text", x = x_loc, y = y_loc, hjust = 0,
             label = paste("# down =", n_downreg))+
    annotate("text", x = x_loc, y = y_loc+up_offset, hjust = 0,
             label = paste("# up =", n_upreg))+
    scale_size_manual(values = c(0.1,1))+
    scale_color_manual(values = c("black", "blue"))+
    geom_point()+
    labs(title = title,
         y = "-log10(FDR)")
}

plot_gene_CPM = function(y, metadata, genes){
  cpm(y) %>%
    as.data.frame() %>%
    mutate(SYMBOL = y$genes$genes) %>%
    filter(SYMBOL %in% genes) %>%
    pivot_longer(-SYMBOL, names_to = "Run", values_to = "CPM") %>%
    left_join(metadata) %>%
    ggplot(aes(x = paste(trt, Run), y = CPM, fill = trt,
               color = trt))+
    # scale_y_log10()+
    geom_bar(stat = "identity")+
    facet_wrap(~SYMBOL, scales = "free_y")+
    theme(axis.text.x = element_blank())+
    labs(x = "sample")
}

plot_cpm_day_v_trt = function(y, gene_list){
  cpm(y) %>%
    as.data.frame() %>%
    mutate(SYMBOL = y$genes$genes) %>%
    filter(SYMBOL %in% gene_list) %>%
    pivot_longer(-SYMBOL, names_to = "Run", values_to = "CPM") %>%
    left_join(metadata) %>%
    filter(!str_detect(trt, "continuous")) %>%
    mutate(SYMBOL = factor(SYMBOL,
                           levels = gene_list)) %>%
    ggplot(aes(x =Day, y = CPM, fill = Treatment ,
               color = Treatment ))+
    geom_point()+
    scale_y_log10()+
    stat_summary(fun = mean, geom = "line", aes(group = Treatment))+
    stat_summary(fun = mean, shape = "-", size =3)+
    facet_wrap(~SYMBOL, scales = "free_y")+
    labs(x = "Day")
}

get_top_25 = function(dgelrt, ranking = "P"){
  # returns a vector of 25 genes in order of max to min. ranking 
  # should be one of c("P", "down", "up")
  if (ranking == "P"){
    ranked_list = as.data.frame(topTags(dgelrt,
                                        n = 25, sort.by = "PValue"))$genes
  }else if (ranking == "down"){
    ranked_list = as.data.frame(dgelrt) %>%
      mutate(FDR = p.adjust(PValue, method = "BH")) %>%
      filter(FDR < 0.05) %>%
      top_n(-25, logFC) %>%
      arrange(logFC) %>%
      pull(genes)
  }else if (ranking == "up"){
    ranked_list = as.data.frame(dgelrt) %>%
      mutate(FDR = p.adjust(PValue, method = "BH")) %>%
      filter(FDR < 0.05) %>%
      top_n(25, logFC) %>%
      arrange(desc(logFC)) %>%
      pull(genes)
  }
  ranked_list
}

get_log_CPM_and_stats_on_geneset = function(dgelist, dgelrt, gene_symbols = NA){
  # This gets the logCPM value for each sample, as well as the statistics, for
  # the symbols in gene_symbols
  # dgelist is commonly called "y" in the vignettes
  # dgelrt is the result of a test like glmQLFTest
  # gene_symbols is a vector of desired genes, in SYMBOL name scheme. If left blank,
  # all are used
  #
  # e.g get_log_CPM_and_stats_on_geneset(y, stats, c("CD226", "CD224"))
  if (is.na(gene_symbols[1])){
    gene_symbols = clusterProfiler::bitr(dgelist$genes_genes,
                                         fromType = "ENSEMBL",
                                         toType = "SYMBOL",
                                         OrgDb = "org.Hs.eg.db")$SYMBOL
  }
  cpm(dgelist, log = TRUE) %>%
    as.data.frame() %>%
    mutate(ENSEMBL = dgelist$genes$genes) %>%
    left_join(clusterProfiler::bitr(.$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL",
                                    OrgDb = "org.Hs.eg.db")) %>%
    left_join(dgelrt %>%
                as.data.frame() %>%
                mutate(FDR = p.adjust(PValue, method = "BH")) %>%
                mutate(ENSEMBL = genes) %>%
                dplyr::select(logFC, FDR, ENSEMBL)) %>%
    filter(SYMBOL %in% gene_symbols) %>%
    pivot_longer(cols = -c(ENSEMBL, SYMBOL, logFC, FDR),
                 names_to = "Run", values_to = "logCPM")
}


myers_metadata_wrangle = function(counts){
  metadata = data.frame(Run = colnames(counts)[7:ncol(counts)]) %>%
    mutate(Day = str_split(Run, "_", simplify = TRUE)[,2],
           Donor = str_split(Run, "_", simplify = TRUE)[,3],
           Treatment = str_split(Run, "_", simplify = TRUE)[,4],
           Gapped = str_split(Run, "_", simplify = TRUE)[,5]) %>%
    mutate(Gapped = case_when(Gapped == "" ~ NA,
                              Gapped == "c" ~ "continuous",
                              Gapped == "g" ~ "gapped")) %>%
    mutate(group = paste(Day, Treatment, sep = ",")) %>%
    mutate(Treatment = factor(Treatment, levels = c("iso", "X"))) %>%
    mutate(Day = factor(Day, levels = c("4", "10"))) %>%
    mutate(trt = paste(Day, Treatment,Gapped, sep = "_")) %>%
    mutate(group_nocon = ifelse(str_detect(Run, "_c"),
                                NA, group))
  metadata
}


myers_load_data_fit_model_no_gapped = function(path_to_counts,
                                               function_for_metadata){
  ## This Loads teh CHURP results, makes the metadata, removes the "gapped"
  ## samples, fits the model, and then returns the counts, the metadata, the 
  ## dgelist, and the fit
  counts = read.table(path_to_counts,
                      header = TRUE)
  metadata = function_for_metadata(counts)
  counts_nondata_cols = names(counts)[1:6]

  
  keep_columns = which(metadata$Treatment_Duration != "gapped")
  metadata = metadata[keep_columns, ]
  counts = counts[, c(1:6, 6 + keep_columns)]
  
  # make edgeR objects
  y = edgeR::DGEList(counts = counts[, 7:ncol(counts)], group = metadata$group,
                     genes = counts$Geneid)

  keepers = edgeR::filterByExpr(y)
  #table(keepers)
  y = y[keepers, , keep.lib.sizes = FALSE]
  y = edgeR::calcNormFactors(y)
  
  # estimate dispersion, so the predFC can use it
  design = model.matrix(~0 + Donor + Day + Treatment, data = metadata)
  
  y = estimateDisp(y, design = design, robust = TRUE)
  
  # make the GLM object
  fit = edgeR::glmQLFit(y, design)
  return(list("metadata" = metadata, 
              "counts" = counts, 
              "dgelist" = y, 
              "fit" = fit))
}

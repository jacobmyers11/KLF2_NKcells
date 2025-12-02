make_ridgeplot = function(ora, dgelrt, xlab = "logFC"){
  # given the result of an ORA test from clusterprofiler and the edgeR F test,
  # makes a ridgeplot. 
  get_ora_data(ora, dgelrt) %>%
    mutate(Description = ifelse(str_length(Description) > 55, 
                                paste(str_sub(Description, 1, 55),"..."),
                                Description)) %>%
    mutate(Description = factor(Description,
                                levels = unique(Description[order(p.adjust,
                                                                  decreasing = TRUE)]))) %>%
    ggplot(aes(x = logFC, y = Description, fill = p.adjust,label = geneID))+
    scale_fill_gradient2(low = "red", mid = "orange", high = "white", midpoint = 0.03)+
    ggridges::geom_density_ridges(jittered_points = TRUE,
                        point_shape = "|",
                        point_size = 2, point_colour = "white",
                        position = ggridges::position_points_jitter(width = 0,
                                                          height = 0,
                                                          yoffset = 0.2))+
    theme(legend.key.height = unit(3, "mm"))+
    labs(x = xlab,
         y = "")
}

make_genesets_per_gene_graph = function(ora, n = 50){
  # this makes a plot of the number of genesets per gene, sorted by the top amount
  # and only showing as many as you want (50 by default)
  # ora is the result of gsea or enricher
  #
  # this also makes the label be the gene sets enriched, this is useful
  # when wrapping the ggplot object with plotly::ggplotly
  #
  # example:
  #  ora = enricher(genes, ...)
  #  graph = make_genesets_per_gene_graph(ora)
  #. plotly::ggplotly(graph)
  ora_df = make_gsea_incidence_matrix(ora) 
  genesets = sapply(1:ncol(ora_df), FUN = function(i){
    paste(rownames(ora_df)[ora_df[,i] > 0], collapse = "\n")
  })
  ora_df = data.frame(gene = colnames(ora_df), 
                      hits = colSums(ora_df),
                      genesets = genesets)
  ora_df %>%
    top_n(n) %>%
    mutate(gene = factor(gene, 
                         levels = gene[order(hits)])) %>%
    ggplot(aes(x = gene, y= hits, label = genesets))+
    geom_bar(stat = "identity")+
    coord_flip()+
    labs(x = "", y = "# gene sets gene is core within")
}


## simplify_go_plus is a hacky way I made to get the benefit of "simplify" while
## using more than just GO terms. It does the enrichment on everything, then
## re-does the enrichment just in GO-BP, simplifies this, and only keeps original
## GO terms if they wre in the simplified version. 
simplify_go_plus = function(genes, y, db){
  # assumes you are going to enrich on GO-BP in addition to other genesets,
  # and you'd like to use "simplify" on just the GO-BP sets. this is messy, but 
  # seems to work. basically, it simplifies a separate GO enrichment, then
  # merges it back with the original enrichment
  
  # the main enrichment
  ora_all = enricher(genes, universe = y$genes$genes, TERM2GENE = db)
  # redo on just GOBP
  ora_go = enrichGO(genes, universe = y$genes$genes, ont = "BP", 
                    OrgDb = org.Hs.eg.db,keyType = "SYMBOL")
  go_simple = simplify(ora_go)
  reduced_go_terms = go_simple %>% as.data.frame %>%
    mutate(term = toupper(Description)) %>%
    mutate(term = paste("GOBP_", term, sep = "")) %>%
    mutate(term = str_replace_all(term, " ", "_")) %>%
    mutate(term = str_replace_all(term, ",", "")) %>%
    mutate(term = str_replace_all(term, "-", "_")) %>%
    pull(term)
  
  go_terms_in_all = ora_all$Description[str_detect(ora_all$Description, "GOBP")]
  other_terms_in_all = ora_all$Description[!str_detect(ora_all$Description, "GOBP")]
  go_keepers = intersect(reduced_go_terms, go_terms_in_all)
  all_keepers = c(go_keepers, other_terms_in_all)
  ora_all %>%
    filter(ID %in% all_keepers)
}

show_categories = function(db){
  # db is the result of msigdbr:msigdbr, which will have gs_cat, gs_subcat, 
  # this shows the unique ones of those
  db %>% 
    dplyr::select(gs_cat, gs_subcat) %>%
    dplyr::group_by(gs_cat, gs_subcat) %>%
    dplyr::filter(dplyr::row_number() == 1) 
}

get_sorted_genelist = function(contrast){
  # contrast is result of glmQLFTest(...)
  # this produces a named vector where the values are logFC
  # and the names are gene names, and returns them sorted in
  # descending order (so higher FC is at the beginning)
  contrast %>%
    as.data.frame() %>%
    dplyr::select(genes, logFC) %>%
    #dplyr::mutate(logFC = -logFC) %>%
    dplyr::arrange(desc(logFC)) %>% 
    dplyr::pull(name = "genes")
}

get_ora_data = function(ora_result, dgelrt){
  # produces a table with enriched terms and their genes along with logFC
  # this is very useful for making a ridgeplot from ora results
  # example:
  #   res = glmQLFTest(...)
  #   ora = enrichGO(...)
  #   dat = get_ora_data(ora, res)
  
  (
    lapply(1:nrow(ora_result), function (x) {
      this_row = as.data.frame(ora_result)[x,]
      genes = str_split(this_row$geneID, "/", simplify= TRUE)[1,]
      this_df = bind_rows(replicate(length(genes), this_row, simplify = FALSE)) %>%
        mutate(genes = genes)
      this_df
    }) %>%
      bind_rows() %>%
      left_join(as.data.frame(dgelrt) %>%
                  dplyr::select(genes, logFC))
  )
}

make_gsea_incidence_matrix = function(gsea_readable){
  # this makes an incidence matrix, where rows are genesets, cols are genes,
  # and the data are 1 if the geneset has the gene in the core, else 0
  if (!(class(gsea_readable)[1] == "data.frame")){
    gsea_readable = as.data.frame(gsea_readable)
  }
  
  # deal with fact that gsea result says "core_enrichment" but
  # ora result says "geneID"
  if ("geneID" %in% colnames(gsea_readable)){
    gsea_readable$core_enrichment = gsea_readable$geneID
  }
  
  genes = sapply(unique(gsea_readable$Description), FUN = function(x){
    stringr::str_split(gsea_readable[gsea_readable$Description == x, "core_enrichment"], "/")
  })
  all_genes = unique(unlist(genes))
  incidence_matrix = matrix(0, nrow = nrow(gsea_readable),
                            ncol = length(all_genes),
                            dimnames = list(gsea_readable$Description,
                                            all_genes))
  # must be a better way than a nested for loop, but it shouldn't be that bad
  # with a few dozen enriched terms
  for (geneset in names(genes)){
    these_genes = genes[[geneset]]
    for (this_gene in these_genes){
      incidence_matrix[geneset, this_gene] = 1
    }
  }
  incidence_matrix
}
reduce_incidence_matrix = function(incidence_matrix, min_colsum = 2){
  # removes columns where colSums are less than min_colsum, e.g. to remove
  # "species" which are only observed once. it then removes any
  # genesets which are no longer represented by any genes
  incidence_matrix = incidence_matrix[ , colSums(incidence_matrix) >= min_colsum]
  incidence_matrix[rowSums(incidence_matrix) > 0, ]
}

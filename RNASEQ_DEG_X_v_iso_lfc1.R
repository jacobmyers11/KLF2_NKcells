####
# This is the main workhorse script for the RNASEQ analysis, which
# uses the objects saved by RNASEQ_create_edgeR_model.R. It calculates
# all the DEG statistics (using FDR < 0.05 and abs(log2FC) > 1), makes plots,
# compares to some other papers' data, and does over-representation analysis.
# It saves figures and some intermediate objects which are used in the 
# report Rmarkdown files. 
#
# To reduce code-reuse, it uses a lot of bespoke figure-making functions i
# the analysis_functions.r file

rm(list = ls())
source("./r_scripts/analysis_functions.r")
load("./intermediate_data/rnaseq_model_min.count_200.Rdat")
library(ggupset)

# contrasts of interest are X vs. iso in day 4, day 10; and change through 
# time of X or iso

# treatment effects within a day
trt_in_day4_contrast = makeContrasts(trt_in_day4_contrast = (trt4_X_NA - trt4_iso_NA),
                                     levels = design)
trt_in_day10_contrast = makeContrasts(trt_in_day10_contrast = (trt10_X_gapped - trt10_iso_gapped),
                                     levels = design)


# time effects within a treatment
X_over_time_contrast = makeContrasts(xot = (trt10_X_gapped - trt4_X_NA),
                                      levels = design)
iso_over_time_contrast = makeContrasts(iot = (trt10_iso_gapped - trt4_iso_NA),
                                     levels = design)

trt_in_day4 = glmQLFTest(fit, contrast = trt_in_day4_contrast)
trt_in_day10 = glmQLFTest(fit, contrast = trt_in_day10_contrast)
x_over_time = glmQLFTest(fit, contrast = X_over_time_contrast)
iso_over_time = glmQLFTest(fit, contrast = iso_over_time_contrast)

# summary tables
get_num_DEGs(trt_in_day4,logFC.thresh = 1) # 623 down, 463 up
get_num_DEGs(trt_in_day10,logFC.thresh = 1)# 208 down, 207 up
get_num_DEGs(x_over_time,logFC.thresh = 1) # 547 down, 487 up
get_num_DEGs(iso_over_time,logFC.thresh = 1) # 163 down, 121 up

data.frame(contrast = c("X vs. iso, day 4",
                        "X vs. iso, day 10",
                        "day 10 vs. day 4, X",
                        "day 10 vs. day 4, iso"),
           up = c(463,207,487,121),
           down = c(623,208,547,163))

# DEGs for upset plot
trt_in_day4_DEGs = get_DEGs(trt_in_day4, logFC.thresh = 1)
trt_in_day10_DEGs = get_DEGs(trt_in_day10, logFC.thresh = 1)
x_over_time_DEGs = get_DEGs(x_over_time, logFC.thresh = 1)
iso_over_time_DEGs = get_DEGs(iso_over_time, logFC.thresh = 1)

DEG_table = data.frame(contrast = "down in X vs. iso on day 4",
                       DEG = trt_in_day4_DEGs$down) %>%
  rbind(data.frame(contrast = "down in X vs. iso on day 10",
                   DEG = trt_in_day10_DEGs$down)) %>%
  rbind(data.frame(contrast = "up in X vs. iso on day 4",
                   DEG = trt_in_day4_DEGs$up)) %>%
  rbind(data.frame(contrast = "up in X vs. iso on day 10",
                   DEG = trt_in_day10_DEGs$up)) %>%
  rbind(data.frame(contrast = "goes up in X over time",
                   DEG = x_over_time_DEGs$up))%>%
  rbind(data.frame(contrast = "goes down in X over time",
                   DEG = x_over_time_DEGs$down))%>%
  rbind(data.frame(contrast = "goes up in iso over time",
                   DEG = iso_over_time_DEGs$up))%>%
  rbind(data.frame(contrast = "goes down in iso over time",
                   DEG = iso_over_time_DEGs$down)) 

(upsetplot = DEG_table %>%
  group_by(DEG) %>%
  summarize(contrasts = list(contrast)) %>%
  ggplot(aes(x=contrasts)) +
  geom_bar() +
  scale_x_upset()+
  labs(y = "# genes", x = ""))
ggsave("./plots/X_vs_iso_lfc1/upsetplot.png",
       plot = upsetplot,
       dpi = 300, width = 7, height = 3)

DEG_table = DEG_table %>%
  mutate(gene = TRUE) %>%
  pivot_wider(names_from = "contrast", values_from = "gene", values_fill = FALSE) 

save(DEG_table, file = "./intermediate_data/X_vs_iso_lfc1/DEG_table.Rdat")

# how are some genes down in X over time, but only that?
# answer: they also go down in iso, and afe often even diff from iso within a time
# point, BUT without a difference that is abs(log2FC) > 1
only_down_in_X = x_over_time_DEGs$down %>%
  setdiff(trt_in_day4_DEGs$up) %>%
  setdiff(trt_in_day4_DEGs$down) %>%
  setdiff(trt_in_day10_DEGs$up) %>%
  setdiff(trt_in_day10_DEGs$down) %>%
  setdiff(x_over_time_DEGs$up) %>%
  setdiff(iso_over_time_DEGs$down) %>%
  setdiff(iso_over_time_DEGs$up)
plot_cpm_day_v_trt(y, only_down_in_X[1:10])+
  facet_wrap(~SYMBOL, scales = "free_y", nrow = 2)# they do go down, and many 
ggsave("./plots/X_vs_iso_lfc1/cpm_xdownovertime.png",
       dpi = 300, width = 10, height = 4)

# Here is the data support for the previous answer: we see negative logFC 
# just smaller than 1, and often significant by FDR
# go down FDR < 0.05, BUT they logFC > -1
iso_over_time %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>%
  filter(genes %in% only_down_in_X[1:10]) %>%
  dplyr::select(genes, logFC, FDR)
# genes      logFC         FDR
# 1     AUNIP -0.7589888 0.089306893
# 2       JUN -0.8166856 0.001040578
# 3     CENPH -0.9314434 0.009336569
# 4   CCDC167 -0.7661740 0.010257054
# 5  HLA-DRB5 -0.1192441 0.918858766
# 6      RPA3 -0.5610599 0.025522496
# 7     DSCC1 -0.8948008 0.035547253
# 8     TIPIN -0.9138099 0.001678934
# 9    TMEM97 -0.8548204 0.001663610
# 10   PMAIP1 -0.8955024 0.001045408

#########
### MD plots, saved for printing in the html
#########
md_trt_in_day4 = (plot_MD(trt_in_day4, up_offset = 0.5, 
                          y_pct = 0.999, abs_lfc_thresh = 1)+
    theme(legend.position = "none")) %>%
    plotly::ggplotly(tooltip = "label")

md_trt_in_day10 = (plot_MD(trt_in_day10, up_offset = 0.5, 
                           y_pct = 0.999, abs_lfc_thresh = 1)+
                     theme(legend.position = "none")) %>%
  plotly::ggplotly(tooltip = "label")

md_X_over_time = (plot_MD(x_over_time, up_offset = 0.5, 
                           y_pct = 0.999, abs_lfc_thresh = 1)+
                     theme(legend.position = "none")) %>%
  plotly::ggplotly(tooltip = "label")

md_iso_over_time = (plot_MD(iso_over_time, up_offset = 0.5, 
                          y_pct = 0.999, abs_lfc_thresh = 1)+
                    theme(legend.position = "none")) %>%
  plotly::ggplotly(tooltip = "label")

save(md_trt_in_day4, 
     md_trt_in_day10,
     md_X_over_time,
     md_iso_over_time,
     file = "./intermediate_data/X_vs_iso_lfc1/md_plots.Rdat")

####
# CPM plots. I will show all the four treatment (day v treat) even when just subsetting
# on a single day's contrast, because that's probably more informative for right now
####

# day 4

top_genes_by_LFCdown = get_top_25(trt_in_day4, "down")
top_genes_by_LFCup = get_top_25(trt_in_day4, "up")

plot_cpm_day_v_trt(y,top_genes_by_LFCdown )
ggsave("./plots/X_vs_iso_lfc1/CPM_trtinday4_top25byLFCdown.png",
       dpi = 300, width = 8, height = 6)

plot_cpm_day_v_trt(y,top_genes_by_LFCup )
ggsave("./plots/X_vs_iso_lfc1/CPM_trtinday4_top25byLFCup.png",
       dpi = 300, width = 8, height = 6)

# day 10 
top_genes_by_LFCdown = get_top_25(trt_in_day10, "down")
top_genes_by_LFCup = get_top_25(trt_in_day10, "up")

plot_cpm_day_v_trt(y,top_genes_by_LFCdown )
ggsave("./plots/X_vs_iso_lfc1/CPM_trtinday10_top25byLFCdown.png",
       dpi = 300, width = 8, height = 6)

plot_cpm_day_v_trt(y,top_genes_by_LFCup )
ggsave("./plots/X_vs_iso_lfc1/CPM_trtinday10_top25byLFCup.png",
       dpi = 300, width = 8, height = 6)

#  effect of time, iso
top_genes_by_LFCdown = get_top_25(iso_over_time, "down")
top_genes_by_LFCup = get_top_25(iso_over_time, "up")

plot_cpm_day_v_trt(y,top_genes_by_LFCdown )
ggsave("./plots/X_vs_iso_lfc1/CPM_isoovertime_top25byLFCdown.png",
       dpi = 300, width = 8, height = 6)

plot_cpm_day_v_trt(y,top_genes_by_LFCup )
ggsave("./plots/X_vs_iso_lfc1/CPM_isoovertime_top25byLFCup.png",
       dpi = 300, width = 8, height = 6)

#  effect of time, X
top_genes_by_LFCdown = get_top_25(x_over_time, "down")
top_genes_by_LFCup = get_top_25(x_over_time, "up")

plot_cpm_day_v_trt(y,top_genes_by_LFCdown )
ggsave("./plots/X_vs_iso_lfc1/CPM_xovertime_top25byLFCdown.png",
       dpi = 300, width = 8, height = 6)

plot_cpm_day_v_trt(y,top_genes_by_LFCup )
ggsave("./plots/X_vs_iso_lfc1/CPM_xovertime_top25byLFCup.png",
       dpi = 300, width = 8, height = 6)

#######
### scatter plots
# first, of log2FC (x - iso) on day 10 vs. day 4
# then, of log2FC(day 10 vs day 4) of X vs. iso
#######

# are the differences between X and iso within a timepoint similar at both timepoints?
day10v4_scatter = trt_in_day10 %>% 
  as.data.frame() %>%
  mutate(day10 = logFC) %>%
  mutate(day10_FDR = p.adjust(PValue, method = "BH")) %>%
  mutate(day10_sig = day10_FDR < 0.05 & abs(day10) > 1) %>%
  dplyr::select(genes, day10, day10_sig) %>%
  left_join(trt_in_day4 %>%
              as.data.frame() %>%
              mutate(day4 = logFC) %>%
              mutate(day4_FDR = p.adjust(PValue, method = "BH")) %>%
              mutate(day4_sig = day4_FDR < 0.05 & abs(day4) > 1) %>%
              dplyr::select(genes, day4, day4_sig)) %>%
  mutate(`significant\ncontrasts:` = ifelse(day10_sig & day4_sig,
                                            "both",
                                            ifelse(day4_sig,
                                                   "day 4 only",
                                                   ifelse(day10_sig, 
                                                          "day 10 only", "neither")))) %>%
  mutate(`significant\ncontrasts:` = factor(`significant\ncontrasts:`,
                                            levels = c("neither", "day 4 only", "day 10 only","both"))) %>%
  ggplot(aes(x = day4, y = day10, label = genes, 
             color = `significant\ncontrasts:`))+
  scale_color_manual(values = c("neither" = "darkgray",
                                "day 4 only" = "blue",
                                "day 10 only" = "green",
                                "both" = "red"))+
  #geom_hline(yintercept = 0, linetype = "dashed")+
  annotate("segment", x = -1, xend = 1, y = -1, yend = -1, linetype = "dashed")+
  annotate("segment", x = -1, xend = 1, y = 1, yend = 1, linetype = "dashed")+
  annotate("segment", x = -1, xend = -1, y = -1, yend = 1, linetype = "dashed")+
  annotate("segment", x = 1, xend = 1, y = -1, yend = 1, linetype = "dashed")+
  #geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(shape = 20, size = 1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlim(-7,10)+
  ylim(-7,10)+
  #geom_text_repel()+
  labs(x = "logFC (X - iso) in day 4",
       y = "logFC (X - iso) in day 10") 
day10v4_scatter
day10v4_scatter = day10v4_scatter %>%
  plotly::ggplotly(tooltip = c("label","x","y"))
day10v4_scatter


# do X and iso have different changes through time?
overtimeXviso_scatter = x_over_time %>% 
  as.data.frame() %>%
  mutate(x = logFC) %>%
  mutate(x_FDR = p.adjust(PValue, method = "BH")) %>%
  mutate(x_sig = x_FDR < 0.05 & abs(x) > 1) %>%
  dplyr::select(genes, x, x_sig) %>%
  left_join(iso_over_time %>%
              as.data.frame() %>%
              mutate(iso = logFC) %>%
              mutate(iso_FDR = p.adjust(PValue, method = "BH")) %>%
              mutate(iso_sig = iso_FDR < 0.05 & abs(iso) > 1) %>%
              dplyr::select(genes, iso, iso_sig)) %>%
  mutate(`significant\ncontrasts:` = ifelse(x_sig & iso_sig,
                                            "both",
                                            ifelse(x_sig,
                                                   "X over time only",
                                                   ifelse(iso_sig, 
                                                          "iso over time only", "neither")))) %>%
  mutate(`significant\ncontrasts:` = factor(`significant\ncontrasts:`,
                                            levels = c("neither", "X over time only", "iso over time only","both"))) %>%
  ggplot(aes(x = iso, y = x, label = genes, 
             color = `significant\ncontrasts:`))+
  scale_color_manual(values = c("neither" = "darkgray",
                                "X over time only" = "blue",
                                "iso over time only" = "green",
                                "both" = "red"))+
  #geom_hline(yintercept = 0, linetype = "dashed")+
  annotate("segment", x = -1, xend = 1, y = -1, yend = -1, linetype = "dashed")+
  annotate("segment", x = -1, xend = 1, y = 1, yend = 1, linetype = "dashed")+
  annotate("segment", x = -1, xend = -1, y = -1, yend = 1, linetype = "dashed")+
  annotate("segment", x = 1, xend = 1, y = -1, yend = 1, linetype = "dashed")+
  #geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(shape = 20, size = 1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlim(-6,6)+
  ylim(-6,6)+
  #geom_text_repel()+
  labs(x = "logFC (day10 - day4) in iso",
       y = "logFC (day10 - day4) in x") 
overtimeXviso_scatter
overtimeXviso_scatter = overtimeXviso_scatter %>%
  plotly::ggplotly(tooltip = c("label","x","y"))
overtimeXviso_scatter
save(day10v4_scatter, overtimeXviso_scatter, 
     file = "./intermediate_data/X_vs_iso_lfc1/scatter_plots.Rdat")

##################
#### comparisons to other datasets
##################

# this sourced script gets together the different datatables from the papers
# supplied to me by Jake
source("./r_scripts/wrangle_other_data.R")

other_scatter = other_data %>%
  left_join(trt_in_day4 %>%
              as.data.frame() %>%
              mutate(x_v_iso_lfc = logFC,
                     x_v_iso_FDR = p.adjust(PValue, method = "BH")) %>%
              dplyr::select(genes, x_v_iso_lfc, x_v_iso_FDR))%>%
  mutate(myers_sig = x_v_iso_FDR < 0.05 & abs(x_v_iso_lfc) > 1,
         other_sig = FDR < 0.05 & abs(logFC) > 1) %>%
  filter(!is.na(myers_sig)) %>%
  mutate(significance = ifelse(myers_sig & other_sig,
                               "both",
                        ifelse(other_sig,
                               "other",
                        ifelse(myers_sig,
                               "Myers", "neither")))) %>%
  ggplot(aes(x = x_v_iso_lfc, y = logFC, 
             color = significance, group = contrast,
             label = genes))+
  geom_point()+
  stat_smooth(method = "lm")+
  facet_wrap(~contrast, scales = "free_y")+
  labs(x = "logFC (X - iso, day 4) (Myers)",
       y = "logFC from faceted dataset")
other_scatter = plotly::ggplotly(other_scatter)
save(other_scatter, file = "./intermediate_data/X_vs_iso_lfc1/other_data_scatter.Rdat")

(upsetplot= other_data %>%
  filter(FDR < 0.05) %>%
  filter(abs(logFC) > 1) %>%
  dplyr::select(genes, contrast) %>%
  rbind(data.frame(contrast = "myers",
                   genes =  c(get_DEGs(trt_in_day4, logFC.thresh = 1)$up,
                              get_DEGs(trt_in_day4, logFC.thresh = 1)$down))) %>%
    group_by(genes) %>%
    summarize(contrasts = list(contrast)) %>%
    ggplot(aes(x=contrasts)) +
    geom_bar() +
    scale_x_upset()+
    labs(y = "# genes", x = ""))
ggsave("./plots/X_vs_iso_lfc1/upsetplot_otherdata.png",
       plot = upsetplot,
       dpi = 300, width = 12, height = 3)

#########
# overrepresentation analysis
#########

library(clusterProfiler)
library(org.Hs.eg.db)
source("./r_scripts/enrichment_functions.r")
# here are the genesets we wish to test:
#genes which were down in X compared to iso at either time 
#genes which were up in X compared to iso at either time
#genes which went down through time uniquely in X 
#genes which went up through time uniquely in X 

# # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009105
# # that paper shows how the number of DEGs affects the number of enriched pathways
# # in a unimodal way; meaning there is a good sweet-spot. I am going to try to find 
# # that sweet spot. 
# results = data.frame()
# thresholds = c(0.25, 0.5, 0.75, seq(from = 1, to = 2, by = 0.05), 2.25, 2.5, 3)
# for (threshold in thresholds){
#   down_in_X = c(get_DEGs(trt_in_day4, logFC.thresh = threshold)$down,
#                 get_DEGs(trt_in_day10, logFC.thresh = threshold)$down) %>%
#     unique()
#   enrichments = enrichGO(down_in_X,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "BP",
#                          universe = y$genes$genes)
#   results = rbind(results,
#                   data.frame(threshold,
#                              n_genes = length(down_in_X),
#                              n_enrichments = nrow(enrichments)))
# }
# a = results %>%
#   ggplot(aes(x = threshold, y = n_enrichments))+
#   geom_line()+
#   geom_point()
# b = results %>%
#   ggplot(aes(x = n_genes, y = n_enrichments))+
#   geom_point()
# ggsave("./plots/X_vs_iso_lfc1/effect_of_lfc_on_ORA_genesets.png",
#        plot = gridExtra::grid.arrange(a,b, nrow = 1),
#        width = 8, height = 2.5)

# well, ok, looks like the threshold of 1 will work. cool!



## save the database downladed in 20230130 using msigdb
load("./intermediate_data/msigdbr_all_homosapiens_20230130.Rdat")

# for exploration, lets do enrichment on the hallmark pathways,
# GO biological processes, and CP Kegg
get_genesets = function(db, string_matches){
  all_descriptions = db$gs_name %>% unique
  genesets = c()
  for (string_to_match in string_matches){
    genesets = c(genesets, 
                 all_descriptions[str_detect(all_descriptions, string_to_match)])
  }
  db %>%
    filter(gs_name %in% genesets)%>%
    dplyr::mutate(term = gs_name,
                  gene = human_gene_symbol) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::select(term, gene)
}
pathway_terms = c("_ERK","_APOPTOSIS","_INFLAMM","_SURVIV","_CYCL")
subdb = get_genesets(humandb, pathway_terms) %>%
  filter(!str_detect(term, "GOBP_")) %>% # to remove redundancy
  filter(!str_detect(term, "HALLMARK_")) %>% # to remove redundancy
  rbind(humandb %>% 
          filter(gs_subcat == "GO:BP" |
                 gs_cat == "H") %>%
          dplyr::mutate(term = gs_name,
                        gene = human_gene_symbol) %>%
          dplyr::filter(!is.na(gene)) %>%
          dplyr::select(term, gene))

## gene lists to test
down_in_X = c(get_DEGs(trt_in_day4, logFC.thresh = 1)$down,
              get_DEGs(trt_in_day10, logFC.thresh = 1)$down) %>%
  unique() # 652 genes
up_in_X = c(get_DEGs(trt_in_day4, logFC.thresh = 1)$up,
            get_DEGs(trt_in_day10, logFC.thresh = 1)$up) %>%
  unique() # 500 genes
down_in_X_overtime = get_DEGs(x_over_time, logFC.thresh = 1)$down %>%
  setdiff(get_DEGs(iso_over_time, logFC.thresh = 1)$down) # 442 genes
up_in_X_overtime = get_DEGs(x_over_time, logFC.thresh = 1)$up %>%
  setdiff(get_DEGs(iso_over_time, logFC.thresh = 1)$up) # 406 genes



## Here we do the ORA. The function simplify_go_plus was made because
# I wanted to take advantage of semantic similarity to simplify the GO term
# enriched, but that doesn't work on the other gene sets, so I do a two-step
# process. It isn't perfect but seems to work

down_ora = simplify_go_plus(down_in_X, y, subdb) # 71
up_ora = simplify_go_plus(up_in_X, y, subdb) # 117
down_in_time_ora = simplify_go_plus(down_in_X_overtime, y, subdb) # 90
up_in_time_ora = simplify_go_plus(up_in_X_overtime, y, subdb) # 16


library(ggridges)

make_ridgeplot(down_ora, trt_in_day4, "log2FC(X - iso, day 4)")
ggsave("./plots/X_vs_iso_lfc1/ridgeplot_downinX.png",
       dpi = 300, width = 8, height = 8* nrow(down_ora) / 60)
make_ridgeplot(up_ora, trt_in_day4, "log2FC(X - iso, day 4)")
ggsave("./plots/X_vs_iso_lfc1/ridgeplot_upinX.png",
       dpi = 300, width = 8, height = 8* nrow(up_ora) / 60)

make_ridgeplot(down_in_time_ora, x_over_time, "log2FC(day10 - day4, X)")
ggsave("./plots/X_vs_iso_lfc1/ridgeplot_downovertimeinX.png",
       dpi = 300, width = 8, height = 8* nrow(down_in_time_ora) / 60)
make_ridgeplot(up_in_time_ora, x_over_time, "log2FC(day10 - day4, X)")
ggsave("./plots/X_vs_iso_lfc1/ridgeplot_upovertimeinX.png",
       dpi = 300, width = 8, height = 8* nrow(up_in_time_ora) / 60)


#### interactive cnetplots -- they end up being too busy, so lets just do top genes
down_plotly = plotly::ggplotly(make_genesets_per_gene_graph(down_ora, 100))
up_plotly = plotly::ggplotly(make_genesets_per_gene_graph(up_ora, 100))
down_in_time_plotly = plotly::ggplotly(make_genesets_per_gene_graph(down_in_time_ora, 100))
up_in_time_plotly = plotly::ggplotly(make_genesets_per_gene_graph(up_in_time_ora, 100))
save(down_plotly, up_plotly, down_in_time_plotly, up_in_time_plotly,
     file = "./intermediate_data/X_vs_iso_lfc1/barplot_genesetspergene_plotly.Rdat")

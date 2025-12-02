####
# Jake's request:
# Second, I was hoping to do GSEA on the Dogra et al lung resident Nk cell data 
# set and the Kumar et al CD8 tissue resident (CD69+) data sets (both comparing 
# to day 10 iso vs X). Depending on how this looks, I might want heatmaps of leading edge genes.

rm(list = ls())
source("./r_scripts/analysis_functions.r")
load("./intermediate_data/rnaseq_model_min.count_200.Rdat")
library(ggupset)
library(ggpmisc) # for easy annotation of R2 and p-val
# contrasts of interest are X vs. iso in day 4, day 10; and change through 
# time of X or iso

# treatment effects within a day
trt_in_day10_contrast = makeContrasts(trt_in_day10_contrast = (trt10_X_gapped - trt10_iso_gapped),
                                     levels = design)
trt_in_day10 = glmQLFTest(fit, contrast = trt_in_day10_contrast)
trt_in_day10_DEGs = get_DEGs(trt_in_day10, logFC.thresh = 1)

# this sourced script gets together the different datatables from the papers
# supplied to me by Jake
source("./r_scripts/wrangle_other_data.R")

library(clusterProfiler)
library(org.Hs.eg.db)
source("./r_scripts/enrichment_functions.r")


##### 
# GSEA of current data vs. Bai et al gene sets
#####
unique(other_data$contrast)
#[1] "Bai_NSCLC_Tex"                  "Bai_HCC_Tex"                    "Bai_CRC_Tex"                    "dogra_lung_v_blood_bright_NK"  
#[5] "dogra_lln_v_blood_bright_NK"    "KUMAR_lung_CD4_cd69p_v_cd69n"   "KUMAR_lung_CD8_cd69p_v_cd69n"   "KUMAR_spleen_CD4_cd69p_v_cd69n"
#[9] "KUMAR_spleen_CD8_cd69p_v_cd69n" "mouse_atlas"                    "tietscher_T_exhaust"    

contrasts_of_interest = c("KUMAR_lung_CD8_cd69p_v_cd69n",
                          "KUMAR_lung_CD8_cd69p_v_cd69n",
                          "KUMAR_spleen_CD8_cd69p_v_cd69n")
term2gene = other_data %>%
  filter(FDR < 0.05) %>%
  filter(contrast %in% contrasts_of_interest) %>%
  mutate(term = ifelse(logFC < 0,
                       paste(contrast, "_DOWN", sep = ""),
                       paste(contrast, "_UP", sep = ""))) %>%
  mutate(gene = genes) %>%
  dplyr::select(term, gene)

day10_sorted = get_sorted_genelist(trt_in_day10)

# day 10
day10_gsea = GSEA(day10_sorted, TERM2GENE = term2gene, eps = 0)

summary(day10_gsea)
# look for significance: (only the significant terms are returned)
as.data.frame(day10_gsea) %>%
  dplyr::select(Description, p.adjust)


# make a "barcode" plot for any of the terms 
term = "KUMAR_lung_CD8_cd69p_v_cd69n_DOWN"
clusterProfiler::gseaplot(day10_gsea,term)
ggsave("./my_barcode_plot_test.pdf",
       dpi = 300, width = 5, height = 8)

# and make a heatmap for that term
leading_edge = day10_gsea %>%
  as.data.frame() %>%
  filter(Description == term) %>%
  pull(core_enrichment) %>%
  str_split(pattern = "/") %>%
  {.[[1]]}

# get the sub-setted matrix of log2-CPM values for just the leading edge genes
cpmy = cpm(y, log = TRUE)
row.names(cpmy) = y$genes$genes
cpmy = cpmy[row.names(cpmy) %in% leading_edge, ]
# and only keep the samples from day 10
day10_samples = metadata$Run[metadata$trt %in% c("10_X_gapped","10_iso_gapped")]
cpmy = cpmy[, day10_samples]
# to make only a few gene names show up, we'll set the rowname to "" to every other gene
genes_of_interest = c("CD47","SOX13","NKG7")
row.names(cpmy)[!(row.names(cpmy) %in% genes_of_interest)] = ""
# make a super-simple d.f. for the annotation (this can be re-used)
annot = metadata %>%
  dplyr::select(Run, Treatment) %>%
  tibble::column_to_rownames("Run")
# make the heatmap
# Note: this clusters the columns, which could result in the treatments 
# not being side-by-side. If you don't want to cluster the columns, add
# cluster_col = FALSE
# here I re-order the matrix to be in treatment order in case cluster_col = FALSE
sample_order = c(colnames(cpmy)[str_detect(colnames(cpmy), "iso")],
                 colnames(cpmy)[str_detect(colnames(cpmy), "X")])
cpmy = cpmy[, sample_order]
(phm = pheatmap::pheatmap(cpmy,
                   annotation_col = annot,
                   show_colnames = FALSE,
                   treeheight_row = 0,
                   #cluster_col = FALSE, # uncomment to not cluster columns
                   border_color = NA)
)
# save
ggsave("./test_phm.pdf",
       plot = phm,
       width = 5, height = 6)

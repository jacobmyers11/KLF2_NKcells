####
# Jake's request:
# First, I was hoping to generate a heatmap for the Bai et al GSEA leading edge 
# genes in the up and down direction for the CRC data set comparison to the Day 10 iso vs X data set. 
# This will be a pretty big heat map, but I’m hoping to annotate just a few genes of interest and leave most genes unannotated. 

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


#########
# overrepresentation analysis
#########

library(clusterProfiler)
library(org.Hs.eg.db)
source("./r_scripts/enrichment_functions.r")


##### 
# GSEA of current data vs. Bai et al gene sets
#####
term2gene = other_data %>%
  filter(FDR < 0.05) %>%
  filter(str_detect(contrast, "Bai_CRC_Tex")) %>%
  mutate(term = ifelse(logFC < 0,
                       paste(contrast, "_DOWN", sep = ""),
                       paste(contrast, "_UP", sep = ""))) %>%
  mutate(gene = genes) %>%
  dplyr::select(term, gene)

day10_sorted = get_sorted_genelist(trt_in_day10)

# day 10
day10_gsea = GSEA(day10_sorted, TERM2GENE = term2gene, eps = 0)
# grab the leading edge genes
term = "Bai_CRC_Tex_UP"
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
genes_of_interest = c("VCAM1","SOX4","CD27")
row.names(cpmy)[!(row.names(cpmy) %in% genes_of_interest)] = ""
# make a super-simple d.f. for the annotation (this can be re-used)
annot = metadata %>%
  dplyr::select(Run, Treatment) %>%
  tibble::column_to_rownames("Run")
# make the heatmap
(phm = pheatmap::pheatmap(cpmy,
                   annotation_col = annot,
                   show_colnames = FALSE,
                   treeheight_row = 0,
                   border_color = NA)
)
# save
ggsave("./test_phm.pdf",
       plot = phm,
       width = 5, height = 6)

# REPEAT, but with z-scores
rowSDs = function(x){
  apply(x, 1, sd)
}
# here is the z-score CPM table ("z scores were calculated for each gene from the log2(CPM) values")
cpm_z = (cpmy - rowMeans(cpmy)) / rowSDs(cpmy)
# make the heatmap
(phm = pheatmap::pheatmap(cpm_z,
                          annotation_col = annot,
                          show_colnames = FALSE,
                          treeheight_row = 0,
                          border_color = NA)
)
# save
ggsave("./test_phm.pdf",
       plot = phm,
       width = 5, height = 6)

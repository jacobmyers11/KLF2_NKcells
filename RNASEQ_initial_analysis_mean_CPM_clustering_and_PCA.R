rm(list = ls())
library(dplyr)
library(tidyr)
library(edgeR)


#########
#
#. This is the "step 0" for RNASeq--looking at the data with PCA and clustering
#. After this script comes RNASEQ_create_edgeR_model.R
#
## This script looks at total reads, mean CPM,
## batch effects (donor), continuous vs. gapped data
#
## Makes some barplots (total reads, mean CPM),
## histogram (CPM), MDS plots (colored by various metadata),
## clustered heatmaps (with and without gapped data)
#
## Saves a table of the number of reads post-filtering
#########
source("./r_scripts/analysis_functions.r")

path_to_counts ="./data/Counts/subread_counts.txt"


####
## 1. Load data and define metadata
####

counts = read.delim("./data/Counts/subread_counts.txt",
                    check.names = FALSE,skip = 1,
                    stringsAsFactors = FALSE,
                    sep = "\t")

# This pulls out the metadata from the fastq file name based on the information
# emailed to me by Jake on 2023-01-09
metadata = data.frame(Run = colnames(counts)[7:ncol(counts)]) %>%
  mutate(Day = str_split(Run, "_", simplify = TRUE)[,2],
         Donor = str_split(Run, "_", simplify = TRUE)[,3],
         Treatment = str_split(Run, "_", simplify = TRUE)[,4],
         Treatment_Duration = str_split(Run, "_", simplify = TRUE)[,5]) %>%
  mutate(Treatment_Duration = case_when(Treatment_Duration == "" ~ "Day 4",
                                        Treatment_Duration == "c" ~ "continuous",
                                        Treatment_Duration == "g" ~ "gapped")) %>%
  mutate(Treatment_Duration = factor(Treatment_Duration,
                                     levels = c("Day 4",
                                                "gapped",
                                                "continuous"))) %>%
  mutate(group = paste(Day, Treatment, sep = ",")) %>%
  mutate(group = factor(group, 
                        levels = c("4,iso", "4,X", "10,iso", "10,X")))

counts_nondata_cols = names(counts)[1:6]


######
# make figure of # of mapped reads per sample, colored by day / trt
######
#
counts %>% 
  pivot_longer(cols = -all_of(counts_nondata_cols), 
               names_to = "Run", values_to = "count") %>%
  left_join(metadata) %>%
  group_by(Run, Day, Donor, Treatment, group, Treatment_Duration) %>%
  summarize(reads = sum(count)) %>%
  ggplot(aes(x = Run, y = reads, fill = group))+
  geom_bar(stat = "identity")+
  scale_fill_discrete(guide = guide_legend(title = "Day, Trt"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x = "sample", y = "total mapped reads",
       title = "mapped reads per sample")
ggsave("./plots/total_reads_by_run_day_and_treatment.png",
       dpi = 300, width = 6, height = 3)

###### Make our edgeR object, which lets us get CPM. 
# check whether our metadata is in the same order as the counts:
all(names(counts)[7:ncol(counts)] == metadata$Run)

#####
# Use edgeR to get normalized CPM values
#####
# note that we'll have to change our group when doing actual tests
y = edgeR::DGEList(counts = counts[, 7:ncol(counts)], group = metadata$group,
                   genes = counts$Geneid)

# how many are lowly expressed?
keepers = edgeR::filterByExpr(y)
table(keepers) # toss 44592 genes, keep 17268
y = y[keepers, , keep.lib.sizes = FALSE]
y = edgeR::calcNormFactors(y)

## Lets look at cpm, to see if expression varies by treatment
cpm(y) %>%
  data.frame() %>%
  mutate(GeneId = y$genes) %>%
  pivot_longer(cols = -GeneId, names_to = "Run", values_to = "cpm") %>%
  left_join(metadata) %>%
  group_by(group, GeneId) %>%
  summarize(cpm = mean(cpm)) %>%
  ggplot(aes(x = cpm))+
  geom_histogram(aes(fill = group),
                 position = position_dodge(),
                 bins = 11)+
  scale_x_log10()+
  labs(x = "CPM", y = "# of genes", title = "mean CPM per gene per treat")
ggsave("./plots/cpm_per_gene_histogram.png",
       dpi = 300, width = 5, height = 3)

cpm(y) %>%
  data.frame() %>%
  mutate(GeneId = y$genes) %>%
  pivot_longer(cols = -GeneId, names_to = "Run", values_to = "cpm") %>%
  left_join(metadata) %>%
  group_by(Run, Day, Treatment, group) %>%
  summarize(expression = mean(cpm)) %>%
  ggplot(aes(x = group, y = expression, shape = Day, color = Treatment))+
  geom_point()+
  scale_color_brewer(type = "qual")+
  stat_summary(fun = "mean", geom = "point", shape = "-", size = 10, color = "black")+
  labs(title = "mean CPM per treat")
ggsave("./plots/mean_cpm_per_run_day_and_trt.png",
       dpi = 300, width = 5, height = 3)

# stats for above plot
cpm(y) %>%
  data.frame() %>%
  mutate(GeneId = y$genes) %>%
  pivot_longer(cols = -GeneId, names_to = "Run", values_to = "cpm") %>%
  left_join(metadata) %>%
  group_by(Run, Day, Treatment, group) %>%
  summarize(expression = mean(cpm)) %>%
  aov(expression ~ Day * Treatment, data = .) %>%
  summary()

####
# Look at principle components plots to determine blocking variables of importance,
# like Donor
####
mds = plotMDS(y, plot = FALSE)
# sanity check, all should be TRUE
rownames(mds$distance.matrix.squared) == metadata$Run
metadata$group == y$samples$group
library(ggrepel)
data.frame(x = mds$x, y = mds$y, 
           group = factor(y$samples$group,
                          levels = c("4,iso", "4,X", "10,iso", "10,X")),
           Run = rownames(y$samples),
           Donor = metadata$Donor) %>%
  ggplot(aes(x = x , y = y, color = Donor, shape = Donor))+
  scale_shape_manual(values = 0:6)+
  geom_point()+
  #geom_text_repel(aes(label = Run), size = 2, color = "gray")+
  labs(x = paste("PC1, ", 100*round(mds$var.explained[1], 4), "%", sep = ""),
       y = paste("PC2, ", 100*round(mds$var.explained[2], 4), "%", sep = ""),
       title = "MDS, color = Donor")
ggsave("./plots/MDS_by_Donor.png",
       dpi = 300, width = 5, height = 3)


data.frame(x = mds$x, y = mds$y, 
           group = factor(y$samples$group,
                          levels = c("4,iso", "4,X", "10,iso", "10,X")),
           Run = rownames(y$samples),
           Duration = metadata$Treatment_Duration) %>%
  ggplot(aes(x = x , y = y, color = Duration, shape = Duration))+
  scale_shape_manual(values = 0:6)+
  geom_point()+
  #geom_text_repel(aes(label = Run), size = 2, color = "gray")+
  labs(x = paste("PC1, ", 100*round(mds$var.explained[1], 4), "%", sep = ""),
       y = paste("PC2, ", 100*round(mds$var.explained[2], 4), "%", sep = ""))
ggsave("./plots/MDS_by_treatment_duration.png",
       dpi = 300, width = 5, height = 3)

data.frame(x = mds$x, y = mds$y, 
           group = factor(y$samples$group,
                          levels = c("4,iso", "4,X", "10,iso", "10,X")),
           Run = rownames(y$samples),
           Duration = metadata$Treatment_Duration) %>%
  ggplot(aes(x = x , y = y, color = group, shape = group))+
  scale_shape_manual(values = 0:6)+
  geom_point(size = 5)+
  #geom_text_repel(aes(label = Run), size = 2, color = "gray")+
  labs(x = paste("PC1, ", 100*round(mds$var.explained[1], 4), "%", sep = ""),
       y = paste("PC2, ", 100*round(mds$var.explained[2], 4), "%", sep = "")) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20, color = "black"),
        axis.title.y=element_text(size=20)) +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.title.x=element_text(size=20)) +
  theme(axis.line=element_line(color = "black"))

ggsave("./plots/MDS_by_treatment_and_day.png",
       dpi = 300, width = 8, height = 6)


######
# Cluster the CPM values and produce a heatmap
######
library(pheatmap)
hm = cpm(y, log = TRUE) %>%
  cor(method = "pearson") %>%
  pheatmap(annotation_col = metadata %>%
             select(Run, group, Treatment_Duration, Donor) %>%
             column_to_rownames("Run") %>%
             mutate(group = factor(group,
                    levels = c( "4,X","4,iso", "10,iso", "10,X"))))
ggsave("./plots/heatmap_and_clustering_of_CPM_correlations.png",
       plot = hm, dpi = 300, width = 9, height = 7)

#####
## Repeat on data with only the gapped condition
#####

keep_columns = which(!str_detect(metadata$Run, "_c"))
metadata_nog = metadata[keep_columns, ]
counts_nog = counts[, c(1:6, 6 + keep_columns)]
y_nog = edgeR::DGEList(counts = counts_nog[, 7:ncol(counts_nog)], 
                   group = metadata_nog$group,
                   genes = counts_nog$Geneid)
keepers = edgeR::filterByExpr(y_nog)
table(keepers) # toss 44895 genes, keep 16965
y_nog = y_nog[keepers, , keep.lib.sizes = FALSE]
y_nog = edgeR::calcNormFactors(y_nog)

mds_y_nog = plotMDS(y_nog, plot = FALSE)
# sanity check, all should be TRUE
rownames(mds_y_nog$distance.matrix.squared) == metadata_nog$Run
metadata_nog$group == y_nog$samples$group
library(ggrepel)
data.frame(x = mds_y_nog$x, y = mds_y_nog$y, 
           group = factor(y_nog$samples$group,
                          levels = c("4,iso", "4,X", "10,iso", "10,X")),
           Run = rownames(y_nog$samples),
           Donor = metadata_nog$Donor) %>%
  ggplot(aes(x = x , y = y, color = group, shape = group))+
  scale_shape_manual(values = 0:6)+
  geom_point()+
  #geom_text_repel(aes(label = Run), size = 2, color = "gray")+
  labs(x = paste("PC1, ", 100*round(mds_y_nog$var.explained[1], 4), "%", sep = ""),
       y = paste("PC2, ", 100*round(mds_y_nog$var.explained[2], 4), "%", sep = ""),
       title = "MDS, color = treatment")
ggsave("./plots/MDS_by_group_no_continuous_day10.png",
       dpi = 300, width = 5, height = 3)

hm = (cpm(y_nog, log = TRUE) %>%
  cor(method = "pearson")) %>%
  pheatmap(annotation_col = metadata_nog %>%
             rownames_to_column() %>%
             select(Run, group, Donor) %>%
             column_to_rownames("Run") %>%
             mutate(group = factor(group,
                                   levels = c( "4,X","4,iso", "10,iso", "10,X"))))
ggsave("./plots/no_continuous_day10_heatmap_and_clustering_of_CPM_correlations.png",
       plot = hm, dpi = 300, width = 9, height = 7)


# NOTE: must do
#    module load /home/lmnp/knut0297/software/modulesfiles/libpng/1.6.34
#    module load libxml2
#    This was all done in R studio server R 4.1 
# install.packages("renv")
# library(renv)
# renv::init(bare = TRUE) # make an empty renv
install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("devtools")
devtools::install_github("GuangchuangYu/ggtree") 
# cran version didn't work with ggplot
BiocManager::install("clusterProfiler")
install.packages("tidyverse")
install.packages("statmod")
install.packages("pheatmap")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("minfi")
BiocManager::install("DMRcate")
BiocManager::install("missMethyl")
install.packages("plotly")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("FlowSorted.Blood.EPIC")



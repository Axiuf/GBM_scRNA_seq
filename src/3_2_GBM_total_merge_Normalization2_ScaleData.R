# Seurat provides an improved method for normalization in the recent update "Seurat v3". Hence we will explore 
# the different preprocessing in two seperate workflows.
# Please see more information in satijalab.org/seurat/.
# Here is the ScaleData workflow.


# Set environment.
source("requirements.R")


# Get Seurat object.
GBM <- readRDS("./data/2_GBM_total_merge_filtered.rds")


# Normalize data & Calculate cell cycle score.
GBM <- NormalizeData(GBM, normalization.method = "LogNormalize", scale.factor = 100000)

load("./data/cycle.rda")
GBM <- CellCycleScoring(GBM, g2m.features = g2m_genes, s.features = s_genes)


# Identification of highly variable features.
GBM <- FindVariableFeatures(GBM, selection.method = "vst", nfeatures = 3000)


# Scaling the data
all.genes <- rownames(GBM)
GBM <- ScaleData(GBM, features = all.genes)


saveRDS(GBM, file = "./data/3_2_GBM_total_merge_filtered_Normalization2_Scale.rds")
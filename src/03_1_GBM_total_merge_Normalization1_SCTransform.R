# Seurat provides an improved method for normalization in the recent update "Seurat v3",
# Hence we will explore the different preprocessing in two seperate workflows.
# Please see more information in satijalab.org/seurat/
# Here is the SCTransform workflow


# Set environment
source("requirements.R")


# Get Seurat object
GBM <- readRDS("./data/02_GBM_total_merge_filtered.rds")


# Calculate cell cycle score
GBM <- NormalizeData(GBM, normalization.method = "LogNormalize", scale.factor = 100000)

load("./data/cycle.rda")
GBM <- CellCycleScoring(GBM, g2m.features = g2m_genes, s.features = s_genes)


# Using SCTransform in Seurat (replace "NormalizeData", "ScaleData", and "FindVariableFeatures")
GBM <- SCTransform(GBM)


saveRDS(GBM, file = "./data/03_1_GBM_total_merge_filtered_Normalization1_SCT.rds")
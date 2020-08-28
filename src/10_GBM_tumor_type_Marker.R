# Set environment.
source("requirements.R")
dir.create("./plots/10_GBM_tumor_type_Marker")
fig_dpi <- 150


# Get Seurat object.
GBM <- readRDS("./data/9_GBM_tumor_type_SCTIntegrate.rds")


# Determine metrics to plot present in GBM@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(GBM, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.1, 
            order = TRUE,
            min.cutoff = 'q15',
            label = TRUE)
ggsave(filename = "metrics_FeaturePlot.tiff", device = "tiff", path = "./plots/10_GBM_tumor_type_Marker", width = 16, height = 21, dpi = fig_dpi)
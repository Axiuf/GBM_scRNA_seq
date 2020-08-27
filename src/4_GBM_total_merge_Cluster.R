# Since SCTransform is a better way to normalize data, we will use the output of 3_1 for the downstream analysis.


# Set environment.
source("requirements.R")
dir.create("./plots/4_GBM_total_merge_Cluster")
fig_dpi <- 150

# Get Seurat object.
# GBM <- readRDS("./data/3_2_GBM_total_merge_filtered_Normalization2_Scale.rds")
GBM <- readRDS("./data/3_1_GBM_total_merge_filtered_Normalization1_SCT.rds")


# Plot variable features.
top10 <- head(VariableFeatures(GBM), 10)
plot1 <- VariableFeaturePlot(GBM, pt.size = 0.5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(filename = "VariableFeaturePlot.tiff", device = "tiff", path = "./plots/4_GBM_total_merge_Cluster", dpi = fig_dpi)


# Perform linear dimensional reduction
GBM <- RunPCA(GBM, npcs = 75)


# Choose the PCs for the downstream analysis
ElbowPlot(GBM, ndims = 75)
ggsave(filename = "PC1to75_ElbowPlot.tiff", device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 12, height = 7, dpi = fig_dpi)

DimPlot(GBM, reduction = "pca", dims = c(1, 2), pt.size = 0.2, group.by = "orig.ident", split.by = "Phase")
ggsave(filename = "PC1&2_DimPlot.tiff", device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 27, height = 7, dpi = fig_dpi)

VizDimLoadings(GBM, dims = 1:4, reduction = "pca")
ggsave(filename = "PC1to4_VizDimLoadings.tiff", device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 12, height = 12, dpi = fig_dpi)

DimHeatmap(GBM, dims = c(1:3, 20:25), cells = 5000, balanced = TRUE, fast = FALSE)
ggsave(filename = "dim1to3&20to25_DimHeatmap.tiff", device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 24, height = 21, dpi = fig_dpi)
DimHeatmap(GBM, dims = c(40:42, 55:60), cells = 5000, balanced = TRUE, fast = FALSE)
ggsave(filename = "dim40to42&55to60_DimHeatmap.tiff", device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 24, height = 21, dpi = fig_dpi)


# Cluster the cells
# Since "resolution" is the key parameter to determine the number of the clusters, we will check a series of values to see the outcome.
# For SCTransform, more PCs is better. However, if you use the 3_2_Scale here, the number of PCs should be set well-founded.
GBM <- FindNeighbors(GBM, dims = 1:60)
cluster_resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2)
GBM <- FindClusters(GBM, resolution = cluster_resolutions)
Idents(GBM) <- "SCT_snn_res.0.2"


# Run non-linear dimensional reduction (UMAP/tSNE)
GBM <- RunUMAP(GBM, dims = 1:60)
GBM <- RunTSNE(GBM, dims = 1:60)


# Plot UMAP/tSNE under different resolutions.
for(cluster_resolution in cluster_resolutions){
  plot1 <- DimPlot(GBM, reduction = "umap", pt.size = 0.1, group.by = paste0("SCT_snn_res.", cluster_resolution), label = TRUE)
  plot2 <- DimPlot(GBM, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = TRUE)
  plot1 + plot2
  ggsave(filename = paste0("umap_res", cluster_resolution, "_DimHeatmap.tiff"), device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 20, height = 7, dpi = fig_dpi)

  plot1 <- DimPlot(GBM, reduction = "tsne", pt.size = 0.1, group.by = paste0("SCT_snn_res.", cluster_resolution), label = TRUE)
  plot2 <- DimPlot(GBM, reduction = "tsne", pt.size = 0.1, group.by = "orig.ident", label = TRUE)
  plot1 + plot2
  ggsave(filename = paste0("tsne_res", cluster_resolution, "_DimHeatmap.tiff"), device = "tiff", path = "./plots/4_GBM_total_merge_Cluster",width = 20, height = 7, dpi = fig_dpi)
}


saveRDS(GBM, file = "./data/4_GBM_total_merge_filtered_SCT_umap&tsne.rds")
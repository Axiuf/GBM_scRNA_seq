# Set environment.
source("requirements.R")
fig_dpi <- 150


# Take args from the bash script.
integrate_dim <- commandArgs(TRUE)
integrate_dim <- as.numeric(integrate_dim)


# Creat output dir.
dir.create("./plots/09_GBM_tumor_type_SCTIntegrate")
plots_dir <- paste0("./plots/09_GBM_tumor_type_SCTIntegrate/", "integrate_dim_", integrate_dim[1])
dir.create(plots_dir)
data_dir <- "./data/09_GBM_tumor_type_SCTIntegrate"
dir.create(data_dir)


# Get the Seurat object and subset malignant cells.
GBM <- readRDS("./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
GBM <- subset(GBM, subset = cell_type == "Malignant cell")


# Rearrange patient information.
GBM$patient <- GBM$orig.ident
GBM$patient[GBM$patient == "P653_jc"] <- "P653"
GBM$patient[GBM$patient == "P673_jc"] <- "P673"

GBM <- subset(GBM, subset = patient != "P887")


# Split object by tumor type and re-normalize.
GBM.list <- SplitObject(GBM, split.by = "tumor_type")

for (i in 1:length(GBM.list)) {
  GBM.list[[i]] <- SCTransform(GBM.list[[i]])
}


# Standerd SCTIntegrate workflow.
GBM.features <- SelectIntegrationFeatures(object.list = GBM.list, nfeatures = 3000)
GBM.list <- PrepSCTIntegration(object.list = GBM.list, anchor.features = GBM.features)


GBM.anchors <- FindIntegrationAnchors(object.list = GBM.list, normalization.method = "SCT", dims = 1:integrate_dim,
                                      anchor.features = GBM.features)
GBM.integrated <- IntegrateData(anchorset = GBM.anchors, normalization.method = "SCT", dims = 1:integrate_dim)
DefaultAssay(GBM.integrated) <- "integrated"


# Standerd cluster workflow
GBM.integrated <- RunPCA(GBM.integrated)

GBM.integrated <- FindNeighbors(GBM.integrated, dims = 1:50)
cluster_resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2)
GBM.integrated <- FindClusters(GBM.integrated, resolution = cluster_resolutions)
Idents(GBM.integrated) <- "integrated_snn_res.0.4"

GBM.integrated <- RunUMAP(GBM.integrated, dims = 1:50)


# Plot UMAP.
for(cluster_resolution in cluster_resolutions){
  plot1 <- DimPlot(GBM.integrated, reduction = "umap", pt.size = 0.1, group.by = paste0("integrated_snn_res.", cluster_resolution), label = TRUE)
  plot2 <- DimPlot(GBM.integrated, reduction = "umap", pt.size = 0.1, group.by = "tumor_type")
  plot1 + plot2
  ggsave(filename = paste0("umap_res", cluster_resolution, "_DimHeatmap.tiff"), device = "tiff", path = plots_dir, width = 16, height = 7, dpi = fig_dpi)
}


# Determine metrics to plot present in GBM@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(GBM.integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.1, 
            order = TRUE,
            min.cutoff = 'q15',
            label = TRUE)
ggsave(filename = "metrics_FeaturePlot.tiff", device = "tiff", path = plots_dir, width = 16, height = 21, dpi = fig_dpi)


saveRDS(GBM.integrated, file = paste0(data_dir, "/", "integrate_dim_", integrate_dim[1], ".rds"))
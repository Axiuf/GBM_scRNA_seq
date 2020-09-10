# First, divide all data into three data sets by their tumor type: "Primary", "Secondary", "Recurrent"
# Then, SCTIntegrate these three data sets


# Set environment
source("requirements.R")


# Take arguments from the bash script
integrate_dim <- commandArgs(TRUE)
integrate_dim <- as.numeric(integrate_dim)


# Create output dir
plots_dir <- "./plots/09_1_GBM_tumor_type_SCTIntegrate"
dir.create(plots_dir)
plots_subdir <- paste0(plots_dir, "/integrate_dim_", integrate_dim[1])
dir.create(plots_subdir)
data_dir <- "./data/09_1_GBM_tumor_type_SCTIntegrate"
dir.create(data_dir)


# Get the Seurat object and subset malignant cells
GBM <- readRDS("./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
GBM <- subset(GBM, subset = cell_type == "Malignant cell")


# Rearrange patient information
GBM$patient <- GBM$orig.ident
GBM$patient[GBM$patient == "P653_jc"] <- "P653"
GBM$patient[GBM$patient == "P673_jc"] <- "P673"

GBM <- subset(GBM, subset = patient != "P887")


# Split object by tumor type and re-normalize
GBM.list <- SplitObject(GBM, split.by = "tumor_type")

for (i in 1:length(GBM.list)) {
  GBM.list[[i]] <- SCTransform(GBM.list[[i]])
}


# Standard SCTIntegrate work flow
GBM.features <- SelectIntegrationFeatures(object.list = GBM.list, nfeatures = 3000)

GBM.list <- PrepSCTIntegration(object.list = GBM.list, anchor.features = GBM.features)

GBM.anchors <- FindIntegrationAnchors(object.list = GBM.list, normalization.method = "SCT", dims = 1:integrate_dim,
                                      anchor.features = GBM.features)
GBM.integrated <- IntegrateData(anchorset = GBM.anchors, normalization.method = "SCT", dims = 1:integrate_dim)
DefaultAssay(GBM.integrated) <- "integrated"


# Standard cluster work flow
GBM.integrated <- RunPCA(GBM.integrated)

GBM.integrated <- FindNeighbors(GBM.integrated, dims = 1:50)
cluster_resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2)
GBM.integrated <- FindClusters(GBM.integrated, resolution = cluster_resolutions)
Idents(GBM.integrated) <- "integrated_snn_res.0.4"

GBM.integrated <- RunUMAP(GBM.integrated, dims = 1:50)


# Plot UMAP
for(cluster_resolution in cluster_resolutions){
  plot1 <- DimPlot(GBM.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", cluster_resolution),
                   pt.size = 0.1, label = TRUE)
  plot2 <- DimPlot(GBM.integrated, reduction = "umap", group.by = "tumor_type", pt.size = 0.1)
  plot3 <- DimPlot(GBM.integrated, reduction = "umap", group.by = "patient", pt.size = 0.1, )
  plot1 + plot2 + plot3
  ggsave(filename = paste0("umap_res", cluster_resolution, "_DimPlot.tiff"),
         device = "tiff", path = plots_subdir, width = 16, height = 7, dpi = fig_dpi)
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
ggsave(filename = "metrics_FeaturePlot.tiff",
       device = "tiff", path = plots_subdir, width = 16, height = 21, dpi = fig_dpi)


saveRDS(GBM.integrated, file = paste0(data_dir, "/", "integrate_dim_", integrate_dim[1], ".rds"))
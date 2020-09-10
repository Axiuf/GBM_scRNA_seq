# When you have got a comprehensive understand of your data, a simplified work flow will help a lot since
# many explorative steps are redundant and only a few key parameters need to be modified.
# Here we add "vars.to.regress = "percent.mt" in data normalization process to regress out the impact of
# mitochondrial genes.
# We also adjust some QC thresholds


# Set environment
source("requirements.R")
plot_dir <- "./plots/06_GBM_total_merge_Simplified_workflow"
dir.create(plot_dir)


# Get Seurat object
GBM <- readRDS("./data/01_Creat_seurat_objects/merged_objects/GBM_total_merge.rds")


# QC
# Add more metrics to meta.data
GBM[["percent.mt"]] <- PercentageFeatureSet(GBM, pattern = "^MT-")
GBM$log10GenesPerUMI <- log10(GBM$nFeature_RNA) / log10(GBM$nCount_RNA)


# Filter cells in the raw data
GBM <- subset(x = GBM, subset = (nCount_RNA >= 1000) & 
                         (nFeature_RNA >= 600) & 
                         (log10GenesPerUMI > 0.80) & 
                         (percent.mt < 10))


# Filter genes in the raw data
counts <- GetAssayData(object = GBM, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
counts <- counts[keep_genes, ]
GBM <- CreateSeuratObject(counts, meta.data = GBM@meta.data)
Idents(GBM) <- "orig.ident"


# Calculate cell cycle score
GBM <- NormalizeData(GBM, normalization.method = "LogNormalize", scale.factor = 100000)

load("./data/cycle.rda")
GBM <- CellCycleScoring(GBM, g2m.features = g2m_genes, s.features = s_genes)


# Using sctransform in Seurat (replace "NormalizeData", "ScaleData", and "FindVariableFeatures")
GBM <- SCTransform(GBM, vars.to.regress = "percent.mt")


# Perform linear dimensional reduction
GBM <- RunPCA(GBM, npcs = 75)


# Cluster the cells
GBM <- FindNeighbors(GBM, dims = 1:60)
cluster_resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2)
GBM <- FindClusters(GBM, resolution = cluster_resolutions)
Idents(GBM) <- "SCT_snn_res.0.2"


# Run non-linear dimensional reduction (UMAP/tSNE)
GBM <- RunUMAP(GBM, dims = 1:60)


# Plot UMAP/tSNE under different resolutions
for(cluster_resolution in cluster_resolutions){
  plot1 <- DimPlot(GBM, reduction = "umap", group.by = paste0("SCT_snn_res.", cluster_resolution),
                   pt.size = 0.1, label = TRUE)
  plot2 <- DimPlot(GBM, reduction = "umap", group.by = "orig.ident",
                   pt.size = 0.1, label = TRUE)
  plot1 + plot2
  ggsave(filename = paste0("umap_res", cluster_resolution, "_DimPlot.tiff"),
         device = "tiff", path = plot_dir, width = 20, height = 7, dpi = fig_dpi)
}


# Determine metrics to plot present in GBM@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(GBM, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.1, 
            order = TRUE,
            min.cutoff = 'q15',
            label = TRUE)
ggsave(filename = "metrics_FeaturePlot.tiff",
       device = "tiff", path = plot_dir, width = 16, height = 21, dpi = fig_dpi)


saveRDS(GBM, file = "./data/06_GBM_total_merge_Simplified_workflow_regress_mt.rds")

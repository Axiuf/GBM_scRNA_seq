# Directly SCTIntegrate all data by patient(non-malignant). Ignore tumor type information in this step
# Using reciprocal PCA here to avoid error in large data set computation

# Set environment
source("requirements.R")


# Take arguments from the bash script
integrate_dim <- commandArgs(TRUE)
integrate_dim <- as.numeric(integrate_dim)


# Create output dir
plots_dir <- "./plots/13_GBM_non_malignant_SCIntegrate"
dir.create(plots_dir)
plots_subdir <- paste0(plots_dir, "/integrate_dim_", integrate_dim[1])
dir.create(plots_subdir)
data_dir <- "./data/13_GBM_non_malignant_SCIntegrate"
dir.create(data_dir)


# Get the Seurat object and subset malignant cells
GBM <- readRDS("./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
GBM <- subset(GBM, subset = cell_type != "Malignant cell")


# Rearrange patient information
GBM$patient <- GBM$orig.ident
GBM$patient[GBM$patient == "P653_jc"] <- "P653"
GBM$patient[GBM$patient == "P673_jc"] <- "P673"

GBM <- subset(GBM, subset = patient != "P887")


# Split object by tumor type and re-normalize
GBM.list <- SplitObject(GBM, split.by = "patient")

for (i in 1:length(GBM.list)) {
  GBM.list[[i]] <- SCTransform(GBM.list[[i]])
}


# Standard SCTIntegrate work flow
GBM.features <- SelectIntegrationFeatures(object.list = GBM.list, nfeatures = 3000)
GBM.list <- PrepSCTIntegration(object.list = GBM.list, anchor.features = GBM.features)

for (i in 1:length(GBM.list)) {
  GBM.list[[i]] <- RunPCA(GBM.list[[i]])
}

GBM.anchors <- FindIntegrationAnchors(object.list = GBM.list, normalization.method = "SCT", dims = 1:integrate_dim,
                                      k.filter = 50, anchor.features = GBM.features, reduction = "rpca")
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
  plot3 <- DimPlot(GBM.integrated, reduction = "umap", group.by = "patient", pt.size = 0.1)
  plot1 + plot2 + plot3
  ggsave(filename = paste0("umap_res", cluster_resolution, "_DimPlot.tiff"),
         device = "tiff", path = plots_subdir, width = 24, height = 7, dpi = fig_dpi)
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


# Find marker
DefaultAssay(GBM) <- "RNA"
Idents(GBM) <- "integrated_snn_res.0.1"

signature_markers <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
                       "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
                       "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Immu_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 12, height = 7, dpi = fig_dpi)

signature_markers <- c("IL7R", "CCR7", "S100A4", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7",
                       "FCER1A", "CST3", "PPBP")
DotPlot(GBM, features = signature_markers, col.min = 0) + RotatedAxis()
ggsave(filename = "Immu2_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 10, height = 7, dpi = fig_dpi)

signature_markers <- c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11", "CD14", "AIF1", "CD3E", "CD3D")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Oligo_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 10, height = 7, dpi = fig_dpi)

# Name new clusters
new.cluster.ids <- c("Mono(microglia)", "Mono(microglia)", "oligodendrocytes", "CD8 T", "Mono(microglia)",
                     "Circulating Mono", "Mono(microglia)", "CD14 Mono", "pDC")
names(new.cluster.ids) <- levels(GBM)
GBM <- RenameIdents(GBM, new.cluster.ids)
GBM$cell_type_annotation <- GBM@active.ident
DimPlot(GBM, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
ggsave(filename = "New_cluster_DimPlot.tiff",
       device = "tiff", path = plots_dir, width = 7, height = 7, dpi = fig_dpi)


# Plot cluster by tumor information
Primary <- as.data.frame(table(GBM$new.cluster.ids[GBM$tumor_type == "Primary"]))
Primary$tumor <- "Primary"
Primary$Freq <- Primary$Freq / sum(Primary$Freq) * 100
Secondary <- as.data.frame(table(GBM$new.cluster.ids[GBM$tumor_type == "Secondary"]))
Secondary$tumor <- "Secondary"
Secondary$Freq <- Secondary$Freq / sum(Secondary$Freq) * 100
Recurrent <- as.data.frame(table(GBM$new.cluster.ids[GBM$tumor_type == "Recurrent"]))
Recurrent$tumor <- "Recurrent"
Recurrent$Freq <- Recurrent$Freq / sum(Recurrent$Freq) * 100

group <- rbind(Primary, Secondary, Recurrent)
names(group) <- c("Cluster_new.cluster.ids", "proportion_by_tumor", "tumor")

group %>% 
  ggplot(aes(x = Cluster_new.cluster.ids, y = proportion_by_tumor, fill = tumor)) + 
  geom_col(alpha = 0.8, position = position_dodge2(), width = 0.8) +
  theme_classic()
ggsave(filename = "tumor_by_cluster0.2_geom_bar.tiff", device = "tiff", path = plots_dir, width = 14, dpi = fig_dpi)

saveRDS(group, paste0(data_dir, "/cluster_by_tumor_geom_bar.rds"))


# Chisq.test of cluster by tumor
Primary <- t(table(GBM$new.cluster.ids[GBM$tumor_type == "Primary"]))
Secondary <- t(table(GBM$new.cluster.ids[GBM$tumor_type == "Secondary"]))
Recurrent <- t(table(GBM$new.cluster.ids[GBM$tumor_type == "Recurrent"]))

group <- rbind(Primary, Recurrent, Secondary)
rownames(group) <- c("Primary", "Recurrent", "Secondary")

chisq.test(group)
E <- chisq.test(group)$expected
O <- chisq.test(group)$observed
(O - E)^2 / E
saveRDS(group, paste0(data_dir, "/cluster_by_tumor_chisq_test.rds"))
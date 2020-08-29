# Set environment.
source("requirements.R")
dir.create("./plots/05_GBM_total_merge_Marker")
fig_dpi <- 150


# Get Seurat object.
GBM <- readRDS("./data/04_GBM_total_merge_filtered_SCT_umap&tsne.rds")


# Add group informations to meta.data, for convenience of marker genes identification.
GBM$cell_type <- NA
GBM$cell_type[GBM$SCT_snn_res.0.2 %in% c(0, 9, 10, 19, 21)] <- "Macrophage / Microglia"
GBM$cell_type[GBM$SCT_snn_res.0.2 == 8] <- "Oligodendrocyte"
GBM$cell_type[GBM$SCT_snn_res.0.2 %in% c(15, 22)] <- "T cell"
GBM$cell_type[!(GBM$SCT_snn_res.0.2 %in% c(0, 9, 10, 19, 21, 8, 15, 22))] <- "Malignant cell"

GBM$tumor_type <- NA
GBM$tumor_type[GBM$orig.ident %in% c("P609", "P673", "P673_jc", "P690", "P828_CD133Neq", "P828_CD133Posi", 
                                     "P843_CD133Neg", "P843_CD133posi", "P873", "P912_left0514", "P912_right0514", "P924")] <- "Primary"
GBM$tumor_type[GBM$orig.ident %in% c("P625", "P653", "P653_jc", "P882", "P887")] <- "Recurrent"
GBM$tumor_type[GBM$orig.ident %in% c("P632", "P645", "P679_yang", "P679_yin", "P711")] <- "Secondary"

GBM$infercnv_type_patient <- GBM$orig.ident
GBM$infercnv_type_patient[GBM$SCT_snn_res.0.2 %in% c(0, 9, 10, 19, 21)] <- "Macrophage / Microglia"
GBM$infercnv_type_patient[GBM$SCT_snn_res.0.2 == 8] <- "Oligodendrocyte"
GBM$infercnv_type_patient[GBM$SCT_snn_res.0.2 %in% c(15, 22)] <- "T cell"

GBM$infercnv_type_tumor_patient <- GBM$infercnv_type_patient
GBM$infercnv_type_tumor_patient[GBM$infercnv_type_tumor_patient %in% c("P609", "P673", "P673_jc", "P690", "P828_CD133Neq", "P828_CD133Posi", 
            "P843_CD133Neg", "P843_CD133posi", "P873", "P912_left0514", "P912_right0514", "P924")] <- paste0("Primary_", 
GBM$infercnv_type_tumor_patient[GBM$infercnv_type_tumor_patient %in% c("P609", "P673", "P673_jc", "P690", "P828_CD133Neq", "P828_CD133Posi", 
            "P843_CD133Neg", "P843_CD133posi", "P873", "P912_left0514", "P912_right0514", "P924")])
GBM$infercnv_type_tumor_patient[GBM$infercnv_type_tumor_patient %in% c("P625", "P653", "P653_jc", "P882", "P887")] <- paste0("Recurrent_",
GBM$infercnv_type_tumor_patient[GBM$infercnv_type_tumor_patient %in% c("P625", "P653", "P653_jc", "P882", "P887")])
GBM$infercnv_type_tumor_patient[GBM$infercnv_type_tumor_patient %in% c("P632", "P645", "P679_yang", "P679_yin", "P711")] <- paste0("Secondary_",
GBM$infercnv_type_tumor_patient[GBM$infercnv_type_tumor_patient %in% c("P632", "P645", "P679_yang", "P679_yin", "P711")])

GBM$infercnv_type_cluster <- GBM$SCT_snn_res.0.2
GBM$infercnv_type_cluster <- paste0("cluster", GBM$infercnv_type_cluster)
GBM$infercnv_type_cluster[GBM$SCT_snn_res.0.2 %in% c(0, 9, 10, 19, 21)] <- "Macrophage / Microglia (cluster 0)"
GBM$infercnv_type_cluster[GBM$SCT_snn_res.0.2 == 8] <- "Oligodendrocyte (cluster 9)"
GBM$infercnv_type_cluster[GBM$SCT_snn_res.0.2 %in% c(15, 22)] <- "T cell (cluster 14)"

metadata <- GBM@meta.data


# https://doi.org/10.1016/j.cell.2019.06.024
# Macrophage
FeaturePlot(GBM, features = c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"), reduction = "umap", min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Macrophage_FeaturePlot.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker",width = 16, height = 21, dpi = fig_dpi)

# T cell
FeaturePlot(GBM, features = c("CD2", "CD3D", "CD3E", "CD3G"), reduction = "umap", min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "T_cell_FeaturePlot.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker",width = 16, height = 14, dpi = fig_dpi)

# Oligodendrocyte
FeaturePlot(GBM, features = c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"), reduction = "umap", min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Oligodendrocyte_FeaturePlot.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker",width = 16, height = 21, dpi = fig_dpi)


# https://doi.org/10.1038/s41467-020-17186-5
# myeloid cell, oligodendrocyte, endothelial cell
FeaturePlot(GBM, features = c("CD34", "ESAM", "CD53", "CD74", "MOG", "MBP"), reduction = "umap", min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Normal_cell_markers_FeaturePlot.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker",width = 16, height = 21, dpi = fig_dpi)


# Determine metrics to plot present in GBM@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(GBM, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.1, 
            order = TRUE,
            min.cutoff = 'q15',
            label = TRUE)
ggsave(filename = "metrics_FeaturePlot.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker", width = 16, height = 21, dpi = fig_dpi)


# Plot ncells by cluster or sample
metadata %>% 
  ggplot(aes(x = SCT_snn_res.0.2, fill = orig.ident)) + 
  geom_bar(position='stack', color='white') +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("nCells")
ggsave(filename = "cluster0.2_by_sample_geom_bar.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker", width = 8, dpi = fig_dpi)

metadata %>% 
  ggplot(aes(x = orig.ident, fill = SCT_snn_res.0.2)) + 
  geom_bar(position='stack', color='white') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("nCells")
ggsave(filename = "sample_by_cluster0.2_geom_bar.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker", width = 8, dpi = fig_dpi)


# Plot by cell type
metadata %>% 
  ggplot(aes(x = orig.ident, fill = cell_type)) + 
  geom_bar(position='stack', color='white') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("nCells")
ggsave(filename = "cell_type_geom_bar.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker", width = 8, dpi = fig_dpi)


# Plot by tumor type
DimPlot(GBM, reduction = "umap", pt.size = 0.1, group.by = "tumor_type", label = TRUE)
ggsave(filename = "tumor_type_DimPlot.tiff", device = "tiff", path = "./plots/05_GBM_total_merge_Marker",width = 8, dpi = fig_dpi)


saveRDS(GBM, "./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")

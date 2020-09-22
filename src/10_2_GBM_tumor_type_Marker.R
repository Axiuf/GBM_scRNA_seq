# Set environment
source("requirements.R")


plots_dir <- "./plots/10_2_GBM_tumor_type_Marker"
dir.create(plots_dir)
data_dir <- "./data/10_2_GBM_tumor_type_Marker"
dir.create(data_dir)


# Get Seurat object
GBM <- readRDS("./data/09_2_GBM_tumor_type_SCTIntegrate/integrate_dim_30.rds")
DefaultAssay(GBM) <- "RNA"
Idents(GBM) <- "integrated_snn_res.0.2"


# Plot Clusters
plot1 <- DimPlot(GBM, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, label = TRUE)
plot2 <- DimPlot(GBM, reduction = "umap", group.by = "tumor_type", pt.size = 0.1)
plot3 <- DimPlot(GBM, reduction = "umap", group.by = "patient", pt.size = 0.1)
plot1 + plot2 + plot3
ggsave(filename = "cluster_tumor_patient_DimPlot.tiff",
       device = "tiff", path = plots_dir, width = 24, height = 7, dpi = fig_dpi)


# Plot metrics
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(GBM, reduction = "umap", 
            features = metrics,
            pt.size = 0.1, 
            order = TRUE,
            min.cutoff = 'q15',
            label = TRUE)
ggsave(filename = "metrics_FeaturePlot.tiff",
       device = "tiff", path = plots_dir, width = 16, height = 21, dpi = fig_dpi)


# Plot percent.mt distribution
metadata <- GBM@meta.data
metadata %>% 
  ggplot(aes(color = tumor_type, x = percent.mt, fill = tumor_type)) + 
  geom_histogram(alpha = 0.3, binwidth = 1, position = "dodge2", color = "white") +
  geom_freqpoly(alpha = 1, binwidth = 1, position = "dodge2") + 
  theme_classic()
ggsave(filename = "percent.mt_geom_freqpoly.tiff",
       device = "tiff", path = plots_dir, width = 10, height = 5, dpi = fig_dpi)

VlnPlot(GBM, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), pt.size = 0.01, ncol = 3)
ggsave(filename = "percent.mt_VlnPlot.tiff",
       device = "tiff", path = plots_dir, width = 15, height = 5, dpi = fig_dpi)


# Plot CD133_patient subset
FeaturePlot(GBM, features = c("PROM1"), reduction = "umap",
            min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "PROM1_FeaturePlot.tiff", device = "tiff", path = plots_dir, dpi = fig_dpi)

GBM_CD133 <- subset(GBM,
                    subset = patient %in% c("P828_CD133Posi", "P828_CD133Neq", "P843_CD133posi", "P843_CD133Neg"))
DimPlot(GBM_CD133, reduction = "umap", split.by = "patient", group.by = "integrated_snn_res.0.2",
        pt.size = 1, label = TRUE)
ggsave(filename = "patient_CD133_DimPlot.tiff", device = "tiff", path = plots_dir, width = 28, dpi = fig_dpi)

Idents(GBM) <- "patient"
VlnPlot(GBM, features = "PROM1", idents = c("P828_CD133Posi", "P828_CD133Neq", "P843_CD133posi", "P843_CD133Neg"),
        pt.size = 1)
ggsave(filename = "patient_CD133_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 7, dpi = fig_dpi)
Idents(GBM) <- "integrated_snn_res.0.2"


# Plot different markers
# https://doi.org/10.1038/s41467-020-17186-5
FeaturePlot(GBM, features = c("TOP2A", "AURKB", "FOXM1", "TYMS", "USP1", "EZH2"), reduction = "umap",
            min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Cell_cycle_markers_FeaturePlot.tiff",
       device = "tiff", path = plots_dir, width = 16, height = 21, dpi = fig_dpi)

signature_markers <- c("APOD", "OLIG2", "STMN1", "DCX", "SOX11", "TNC", "CD44", "S100A10",
                       "VIM", "HLA-A", "APOE", "HSPA1B", "DNAJB1", "HSPA6")
FeaturePlot(GBM, features = signature_markers, reduction = "umap",
            min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Cell_type_markers_FeaturePlot.tiff",
       device = "tiff", path = plots_dir, width = 30, height = 30, dpi = fig_dpi)


# Plot cluster by tumor information
Primary <- as.data.frame(table(GBM$integrated_snn_res.0.2[GBM$tumor_type == "Primary"]))
Primary$tumor <- "Primary"
Primary$Freq <- Primary$Freq / sum(Primary$Freq) * 100
Secondary <- as.data.frame(table(GBM$integrated_snn_res.0.2[GBM$tumor_type == "Secondary"]))
Secondary$tumor <- "Secondary"
Secondary$Freq <- Secondary$Freq / sum(Secondary$Freq) * 100
Recurrent <- as.data.frame(table(GBM$integrated_snn_res.0.2[GBM$tumor_type == "Recurrent"]))
Recurrent$tumor <- "Recurrent"
Recurrent$Freq <- Recurrent$Freq / sum(Recurrent$Freq) * 100

group <- rbind(Primary, Secondary, Recurrent)
names(group) <- c("Cluster_integrated_snn_res.0.2", "proportion_by_tumor", "tumor")

group %>% 
  ggplot(aes(x = Cluster_integrated_snn_res.0.2, y = proportion_by_tumor, fill = tumor)) + 
  geom_col(alpha = 0.8, position = position_dodge2(), width = 0.8) +
  theme_classic()
ggsave(filename = "tumor_by_cluster0.2_geom_bar.tiff", device = "tiff", path = plots_dir, width = 14, dpi = fig_dpi)

saveRDS(group, paste0(data_dir, "/cluster_by_tumor_geom_bar.rds"))

# Chisq.test of cluster by tumor
Primary <- t(table(GBM$integrated_snn_res.0.2[GBM$tumor_type == "Primary"]))
Secondary <- t(table(GBM$integrated_snn_res.0.2[GBM$tumor_type == "Secondary"]))
Recurrent <- t(table(GBM$integrated_snn_res.0.2[GBM$tumor_type == "Recurrent"]))

group <- rbind(Primary, Recurrent, Secondary)
rownames(group) <- c("Primary", "Recurrent", "Secondary")

chisq.test(group)
E <- chisq.test(group)$expected
O <- chisq.test(group)$observed
(O - E)^2 / E
saveRDS(group, paste0(data_dir, "/cluster_by_tumor_chisq_test.rds"))


# Identify cluster 9&10
GBM_total_merge <- readRDS("./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
GBM_total_integrate <- readRDS("./data/09_3_GBM_tumor_type_SCTIntegrate/integrate_dim_30.rds")

cluster9 <- rownames(GBM@meta.data)[GBM@meta.data$integrated_snn_res.0.2 == 9]
cluster10 <- rownames(GBM@meta.data)[GBM@meta.data$integrated_snn_res.0.2 == 10]

plot_1 <- DimPlot(GBM_total_merge, reduction = "umap", group.by = "SCT_snn_res.0.01", pt.size = 0.1, label = TRUE)
plot_2 <- DimPlot(GBM_total_merge, reduction = "umap", cells.highlight = cluster9,
                  pt.size = 0.1, sizes.highlight = 0.2) + NoLegend()
plot_3 <- DimPlot(GBM_total_merge, reduction = "umap", cells.highlight = cluster10,
                  pt.size = 0.1, sizes.highlight = 0.2) + NoLegend()
plot_1 + plot_2 + plot_3
ggsave(filename = "GBM_total_merge_tumor_cluster9&10_DimPlot.tiff", 
       device = "tiff", path = plots_dir, width = 21, dpi = fig_dpi)


plot_1 <- DimPlot(GBM_total_integrate, reduction = "umap", group.by = "integrated_snn_res.0.1",
                  pt.size = 0.1, label = TRUE)
plot_2 <- DimPlot(GBM_total_integrate, reduction = "umap", cells.highlight = cluster9,
                  pt.size = 0.1, sizes.highlight = 0.2) + NoLegend()
plot_3 <- DimPlot(GBM_total_integrate, reduction = "umap", cells.highlight = cluster10,
                  pt.size = 0.1, sizes.highlight = 0.2) + NoLegend()
plot_1 + plot_2 + plot_3
ggsave(filename = "GBM_total_integrate_tumor_cluster9&10_DimPlot.tiff", 
       device = "tiff", path = plots_dir, width = 21, dpi = fig_dpi)


plot_1 <- DimPlot(GBM, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, label = TRUE)
plot_2 <- DimPlot(GBM, reduction = "umap", cells.highlight = cluster9,
                  pt.size = 0.1, sizes.highlight = 0.2) + NoLegend()
plot_3 <- DimPlot(GBM, reduction = "umap", cells.highlight = cluster10,
                  pt.size = 0.1, sizes.highlight = 0.2) + NoLegend()
plot_1 + plot_2 + plot_3
ggsave(filename = "GBM_tumor_cluster9&10_DimPlot.tiff", 
       device = "tiff", path = plots_dir, width = 21, dpi = fig_dpi)


# Find all markers
all_markers <- FindAllMarkers(GBM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
saveRDS(all_markers, paste0(data_dir, "/all_markers.rds"))

all_markers_name <- rownames(all_markers)
GBM <- ScaleData(GBM, features = all_markers_name)
saveRDS(GBM, paste0(data_dir, "/10_2_GBM_tumor_type_Marker.rds"))

top20 <- all_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(GBM, features = top20$gene)
ggsave(filename = "all_markers_DoHeatmap.pdf",
       device = "pdf", path = plots_dir, width = 24, height = 24, dpi = fig_dpi)


# Plot cluster markers
DotPlot(GBM, features = "MKI67") + RotatedAxis()
ggsave(filename = "Progenitor_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 4, height = 7, dpi = fig_dpi)

signature_markers <- c("OLIG2", "OLIG1", "APOD", "PLLP", "PLP1", "OMG")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Oligo_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 8, height = 7, dpi = fig_dpi)

signature_markers <- c("GFAP", "S100A10", "HOPX", "CD44")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Astro_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 6, height = 7, dpi = fig_dpi)

signature_markers <- c("STMN1", "SOX11", "DCX", "NSG1")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Neuro_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 6, height = 7, dpi = fig_dpi)

signature_markers <- c("GAPDH", "VEGFA", "HILPDA", "ADM")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Hypoxia_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 6, height = 7, dpi = fig_dpi)

signature_markers <- c("AIF1", "FCGR3A", "MS4A7", "CD68")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Microglia_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 6, height = 7, dpi = fig_dpi)

signature_markers <- c("ISG15", "IFI6", "MX1", "IFIT1")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Interferon_activated_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 6, height = 7, dpi = fig_dpi)

signature_markers <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
                       "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
                       "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Immu_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 12, height = 7, dpi = fig_dpi)

signature_markers <- c("HSPA1B", "DNAJB1", "HSPB1")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Heat_responsed_activated_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 5, height = 7, dpi = fig_dpi)

signature_markers <- c("HLA-A", "APOE", "BCL3", "VIM", "CD44", "ANXA1", "ITGB1", "TGFBI", "LOX", "COL1A2", "VDR",
                       "IL6", "MMP7", "ANXA2", "CHI3L1", "RELB", "TLR2", "TLR4", "CASP1")
DotPlot(GBM, features = signature_markers) + RotatedAxis()
ggsave(filename = "Mes_DotPlot.tiff",
       device = "tiff", path = plots_dir, width = 12, height = 7, dpi = fig_dpi)


# Name new clusters
new.cluster.ids <- c("Oligo", "Astrocytic", "Low quality cell", "Progenitor(S)",
                     "Hypoxia", "Neuronal", "Progenitor(G2/M)",
                     "Low quality cell", "Heat responsed", "Microglia(malignant)", "Interferon activated")
names(new.cluster.ids) <- levels(GBM)
GBM <- RenameIdents(GBM, new.cluster.ids)
DimPlot(GBM, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
ggsave(filename = "New_cluster_DimPlot.tiff",
       device = "tiff", path = plots_dir, width = 7, height = 7, dpi = fig_dpi)
# Set environment
source("requirements.R")

plots_dir <- "./plots/10_1_GBM_tumor_type_Marker"
dir.create(plots_dir)


# Get Seurat object
GBM <- readRDS("./data/09_1_GBM_tumor_type_SCTIntegrate/integrate_dim_30.rds")
Idents(GBM) <- "integrated_snn_res.0.2"


# Plot Clusters
plot1 <- DimPlot(GBM, reduction = "umap", group.by = "integrated_snn_res.0.2", pt.size = 0.1, label = TRUE)
plot2 <- DimPlot(GBM, reduction = "umap", group.by = "tumor_type", pt.size = 0.1)
plot3 <- DimPlot(GBM, reduction = "umap", group.by = "patient", pt.size = 0.1)
plot1 + plot2 + plot3
ggsave(filename = "cluster_tumor_patient_dimplot.tiff",
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
ggsave(filename = "tumor_by_cluster0.2_geom_bar.tiff",
       device = "tiff", path = plots_dir, width = 14, dpi = fig_dpi)


DefaultAssay(GBM) <- "RNA"
# Plot different markers
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
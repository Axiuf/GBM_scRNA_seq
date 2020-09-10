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
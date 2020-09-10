# Set environment
source("requirements.R")
library(readxl)
library(clusterProfiler)
library(DOSE)
library(enrichplot)


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


# Plot different markers.
# https://doi.org/10.1038/s41467-020-17186-5
FeaturePlot(GBM, features = c("TOP2A", "AURKB", "FOXM1", "TYMS", "USP1", "EZH2"), reduction = "umap",
            min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Cell_cycle_markers_FeaturePlot.tiff",
       device = "tiff", path = plots_dir,width = 16, height = 21, dpi = fig_dpi)

signature_markers <- c("APOD", "OLIG2", "STMN1", "DCX", "SOX11", "TNC", "CD44", "S100A10",
                       "VIM", "HLA-A", "APOE", "HSPA1B", "DNAJB1", "HSPA6")
FeaturePlot(GBM, features = signature_markers, reduction = "umap",
            min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "Cell_type_markers_FeaturePlot.tiff",
       device = "tiff", path = plots_dir,width = 30, height = 30, dpi = fig_dpi)


# Plot meta module score by PercentageFeatureSet
# https://doi.org/10.1016/j.cell.2019.06.024
meta_module <- read_excel("data/1-s2.0-S0092867419306877-mmc2.xlsx", skip = 4)
for (i in 1:8){
  score_name <- paste0("meta_module_", names(meta_module)[i])
  module_features <- meta_module[[i]]
  module_features <- module_features[module_features %in% row.names(GBM)]
  print(length(module_features))
  GBM[[score_name]] <- PercentageFeatureSet(GBM, features = module_features)
}
GBM[["meta_module_MES"]] <- GBM[["meta_module_MES2"]] + GBM[["meta_module_MES1"]]
GBM[["meta_module_NPC"]] <- GBM[["meta_module_NPC2"]] + GBM[["meta_module_NPC1"]]

module_names <- c("meta_module_MES", "meta_module_AC", "meta_module_OPC", 
                  "meta_module_NPC", "meta_module_G1.S", "meta_module_G2.M")
FeaturePlot(GBM, features = module_names, reduction = "umap",
            min.cutoff = 'q15', pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "meta_module_FeaturePlot.tiff",
       device = "tiff", path = plots_dir,width = 14, height = 21, dpi = fig_dpi)

VlnPlot(GBM, features = module_names, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_VlnPlot.tiff",
       device = "tiff", path = plots_dir, width = 21, height = 21, dpi = fig_dpi)

VlnPlot(GBM, features = "meta_module_MES", y.max = 20, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_MES_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "meta_module_AC", y.max = 15, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_AC_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "meta_module_OPC", y.max = 7.5, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_OPC_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "meta_module_NPC", y.max = 15, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_NPC_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "meta_module_G1.S", y.max = 5, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_G1.S_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "meta_module_G2.M", y.max = 5, pt.size = 0.01) + NoLegend()
ggsave(filename = "meta_module_G2.M_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)


# Plot TCGA score
# http://dx.doi.org/10.1016/j.ccell.2017.06.003
mmc2 <- read_excel("data/mmc2.xlsx", sheet = "Subtype Genes")
for (i in 1:3){
  score_name <- paste0("TCGA_", names(mmc2)[i])
  TCGA_features <- mmc2[[i]]
  TCGA_features <- TCGA_features[TCGA_features %in% row.names(GBM)]
  print(length(TCGA_features))
  GBM[[score_name]] <- PercentageFeatureSet(GBM, features = TCGA_features)
}

TCGA_names <- c("TCGA_Mesenchymal", "TCGA_Proneural", "TCGA_Classical")

FeaturePlot(GBM, features = TCGA_names, reduction = "umap",
            min.cutoff = 'q15', ncol = 3, pt.size = 0.1, label = TRUE, order = TRUE)
ggsave(filename = "TCGA_FeaturePlot.tiff", device = "tiff", path = plots_dir,width = 21, height = 7, dpi = fig_dpi)

VlnPlot(GBM, features = TCGA_names, pt.size = 0.01) + NoLegend()
ggsave(filename = "TCGA_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 15, height = 7, dpi = fig_dpi)

VlnPlot(GBM, features = "TCGA_Mesenchymal", y.max = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = "TCGA_Mesenchymal_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "TCGA_Proneural", y.max = 1.5, pt.size = 0.01) + NoLegend()
ggsave(filename = "TCGA_Proneural_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)
VlnPlot(GBM, features = "TCGA_Classical", y.max = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = "TCGA_Classical_VlnPlot.tiff", device = "tiff", path = plots_dir, width = 5, dpi = fig_dpi)


# Plot cluster by tumor information.
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


# Chisq.test of cluster by tumor.
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


# Find markers
cluster2.markers <- FindMarkers(GBM, ident.1 = 2, min.pct = 0.25)
saveRDS(cluster2.markers, paste0(data_dir, "/cluster2.markers.rds"))

cluster0.markers <- FindMarkers(GBM, ident.1 = 0, min.pct = 0.25)
saveRDS(cluster0.markers, paste0(data_dir, "/cluster0.markers.rds"))

cluster1.markers <- FindMarkers(GBM, ident.1 = 1, min.pct = 0.25)
saveRDS(cluster1.markers, paste0(data_dir, "/cluster1.markers.rds"))

cluster5.markers <- FindMarkers(GBM, ident.1 = 5, min.pct = 0.25)
saveRDS(cluster5.markers, paste0(data_dir, "/cluster5.markers.rds"))

cluster7.markers <- FindMarkers(GBM, ident.1 = 7, min.pct = 0.25)
saveRDS(cluster7.markers, paste0(data_dir, "/cluster7.markers.rds"))


cluster_markers <- cluster7.markers

cluster_markers <- dplyr::arrange(cluster_markers["avg_logFC"], -avg_logFC)
cluster_markers <- rownames(cluster_markers)[abs(cluster_markers$avg_logFC) > log(2)]
gene.df <- bitr(cluster_markers, fromType = "SYMBOL", 
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
cluster_markers <- gene.df$ENTREZID
cluster_markers <- cluster_markers[!duplicated(cluster_markers)]


ego1 <- enrichGO(gene          = cluster_markers,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

ego2 <- enrichGO(gene          = cluster_markers,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)

p1 <- dotplot(ego1, showCategory=30)
p2 <- dotplot(ego2, showCategory=30)
p1 + p2
ggsave(filename = "cluster7_GO.tiff", device = "tiff", path = plots_dir, width = 16, height = 7, dpi = fig_dpi)

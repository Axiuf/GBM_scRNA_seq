# Set environment
source("requirements.R")
library(readxl)

plots_dir <- "./plots/11_GBM_tumor_type_Signature_score"
dir.create(plots_dir)


# Get Seurat object
GBM <- readRDS("./data/09_2_GBM_tumor_type_SCTIntegrate/integrate_dim_30.rds")
DefaultAssay(GBM) <- "RNA"
Idents(GBM) <- "integrated_snn_res.0.2"


# Plot meta module score
# https://doi.org/10.1016/j.cell.2019.06.024
meta_module <- read_excel("data/supplementaries/1-s2.0-S0092867419306877-mmc2.xlsx", skip = 4)
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
mmc2 <- read_excel("data/supplementaries/mmc2.xlsx", sheet = "Subtype Genes")
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


saveRDS(GBM, "./data/11_GBM_tumor_type_Signature_score.rds")
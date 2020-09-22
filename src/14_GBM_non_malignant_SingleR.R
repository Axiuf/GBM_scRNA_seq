# Set environment
source("requirements.R")
library(SingleR)


# Create output dir
plots_dir <- "./plots/14_GBM_non_malignant_SingleR"
dir.create(plots_dir)
data_dir <- "./data/14_GBM_non_malignant_SingleR"
dir.create(data_dir)


# Get Seurat object
GBM <- readRDS("./data/13_GBM_non_malignant_SCIntegrate/integrate_dim_30.rds")
DefaultAssay(GBM) <- "RNA"
Idents(GBM) <- "integrated_snn_res.0.1"


# Prepare test data and reference data
GBM_sce <- as.SingleCellExperiment(GBM)
MoIm <- MonacoImmuneData()


# SingleR
GBM_pred <- SingleR(test = GBM_sce, ref = MoIm, labels = MoIm$label.main)

saveRDS(GBM_pred, paste0(data_dir, "/GBM_pred.rds"))
saveRDS(GBM@meta.data, paste0(data_dir, "/meta.rds"))


plotScoreHeatmap(GBM_pred)

GBM[["SingleR.labels"]] <- GBM_pred$labels

plot1 <- DimPlot(GBM, reduction = "umap", group.by = "integrated_snn_res.0.1",
                 pt.size = 0.1, label = TRUE)
plot2 <- DimPlot(GBM, reduction = "umap", group.by = "SingleR.labels", pt.size = 0.1, label = TRUE)
plot1 + plot2
ggsave(filename = "SingleR.tiff",
       device = "tiff", path = plots_dir, width = 16, height = 7, dpi = fig_dpi)
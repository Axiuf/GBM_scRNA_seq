# Set environment
source("requirements.R")
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)


plots_dir <- "./plots/12_GBM_tumor_type_GO_enrichment"
dir.create(plots_dir)
data_dir <- "./data/12_GBM_tumor_type_GO_enrichment"
dir.create(data_dir)


# Get Seurat object
GBM <- readRDS("./data/11_GBM_tumor_type_Signature_score.rds")
DefaultAssay(GBM) <- "RNA"
Idents(GBM) <- "integrated_snn_res.0.2"


# Find markers
for(i in 0:10){
  markers <- FindMarkers(GBM, ident.1 = i, min.pct = 0.25)
  saveRDS(markers, paste0(data_dir, "/cluster", i, "_markers.rds"))
  
  
  markers <- dplyr::arrange(markers["avg_logFC"], -avg_logFC)
  markers <- rownames(markers)[markers$avg_logFC > log(2)]
  markers_IDs <- bitr(markers, fromType = "SYMBOL", 
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  markers <- markers_IDs$ENTREZID
  markers <- markers[!duplicated(markers)]
  
  
  ego1 <- enrichGO(gene          = markers,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
  
  ego2 <- enrichGO(gene          = markers,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
  
  p1 <- dotplot(ego1, showCategory = 30)
  p2 <- dotplot(ego2, showCategory = 30)
  p1 + p2
  ggsave(filename = paste0("cluster", i, "_GO.tiff"),
         device = "tiff", path = plots_dir, width = 16, height = 7, dpi = fig_dpi)
}
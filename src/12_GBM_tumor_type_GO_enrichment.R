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


# Find universe gene
universe_gene <- rownames(GBM)
saveRDS(universe_gene, paste0(data_dir, "/universe_gene.rds"))

universe_gene <- readRDS(paste0(data_dir, "/universe_gene.rds"))
universe_IDs <- bitr(universe_gene, fromType = "SYMBOL", 
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)
universe_gene <- universe_IDs$ENTREZID


# Find markers
for(i in 0:10){
  markers <- FindMarkers(GBM, ident.1 = i, min.pct = 0.25)
  saveRDS(markers, paste0(data_dir, "/cluster", i, "_markers.rds"))
}

markers <- FindMarkers(GBM, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
saveRDS(markers, paste0(data_dir, "/cluster0to1_markers.rds"))
markers <- FindMarkers(GBM, ident.1 = 1, ident.2 = 0, min.pct = 0.25)
saveRDS(markers, paste0(data_dir, "/cluster1to0_markers.rds"))


# Go enrichment
for(i in 0:10){
  # Prepare for input
  markers <- readRDS(paste0(data_dir, "/cluster", i, "_markers.rds"))
  markers <- dplyr::arrange(markers, -avg_logFC)
  
  markers_name <- rownames(markers)
  markers_IDs <- bitr(markers_name, fromType = "SYMBOL", 
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  markers_name <- markers_IDs$SYMBOL
  
  markers_List <- markers$avg_logFC[rownames(markers) %in% markers_name]
  names(markers_List) <- markers_IDs$ENTREZID
  markers_name <- names(markers_List)[markers_List > 0]
  
  # GO ORA
  ego1 <- enrichGO(gene          = markers_name,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
  ego1 <- simplify(ego1)
  dotplot(ego1)
  ggsave(filename = paste0("cluster", i, "_enrichGO_dotplot.tiff"),
         device = "tiff", path = plots_dir, width = 8, height = 6, dpi = fig_dpi)
  emapplot(ego1, showCategory = 30, layout = "kk")
  ggsave(filename = paste0("cluster", i, "_enrichGO_emapplot.tiff"),
         device = "tiff", path = plots_dir, width = 8, height = 8, dpi = fig_dpi)
  
  # GO GSEA
  ego2 <- gseGO(geneList     = markers_List,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                pvalueCutoff = 0.05)
  ego2 <- simplify(ego2)
  emapplot(ego2, showCategory = 30, layout = "kk", color = "NES")
  ggsave(filename = paste0("cluster", i, "_gseGO_emapplot.tiff"),
         device = "tiff", path = plots_dir, width = 8, height = 8, dpi = fig_dpi)
  ridgeplot(ego2, showCategory = 20)
  ggsave(filename = paste0("cluster", i, "_gseGO_ridgeplot.tiff"),
         device = "tiff", path = plots_dir, width = 12, height = 10, dpi = fig_dpi)
  
  # KEGG ORA
  kk <- gseKEGG(geneList      = markers_List,
                 minGSSize    = 1,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
  
  ridgeplot(kk, showCategory = 20)
  ggsave(filename = paste0("cluster", i, "_gseKEGG_ridgeplot.tiff"),
         device = "tiff", path = plots_dir, width = 12, height = 10, dpi = fig_dpi)
}
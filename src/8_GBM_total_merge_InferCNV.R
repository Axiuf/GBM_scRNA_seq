# Set environment.
library(Seurat)
library(infercnv)
dir.create("./data/8_GBM_total_merge_InferCNV")


# Get the Seurat object.
GBM <- readRDS("./data/5_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
Idents(GBM) <- "infercnv_type_tumor_patient"


# Rearrange the order of groups.
GBM.list <- SplitObject(GBM, split.by = "infercnv_type_tumor_patient")
GBM.list <- GBM.list[sort(names(GBM.list))]
GBM <- merge(GBM.list[[1]], GBM.list[c(2:25)])


# Get counts matrix.
counts_matrix <- GetAssayData(GBM, slot = "counts")


# Get cell annotations.
cell_Annotations <- GBM@meta.data["infercnv_type_tumor_patient"]
write.table(cell_Annotations, file = "./data/8_GBM_total_merge_InferCNV/cell_Annotations.txt", quote = F, sep= "\t", col.names = FALSE)


# InferCNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = "./data/8_GBM_total_merge_InferCNV/cell_Annotations.txt",
                                    gene_order_file = "./data/gencode_v19_gen_pos.complete.txt",
                                    delim = "\t",
                                    ref_group_names = c("Macrophage / Microglia", "Oligodendrocyte", "T cell"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = "./data/8_GBM_total_merge_InferCNV", 
                             window_length = 101,
                             cluster_by_groups = TRUE, 
                             denoise = TRUE,
                             HMM = TRUE, 
                             save_rds = TRUE)
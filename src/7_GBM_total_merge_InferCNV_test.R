# Considering that this is a explorative step to assess the heterogeneity within tumors and different 
# clusters, I provide a test script for convenience of customization and parallelization.


# Set environment.
library(Seurat)
library(infercnv)


# Take args from the bash script.
win_len <- commandArgs(TRUE)
win_len <- as.numeric(win_len)


# Creat output dir.
dir.create("./data/7_GBM_total_merge_InferCNV_test")
subdir_path <- paste0("./data/7_GBM_total_merge_InferCNV_test/", "window_length_", win_len[1])
dir.create(subdir_path)


# Get the Seurat object and downsample.
GBM <- readRDS("./data/5_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
Idents(GBM) <- "infercnv_type_tumor_patient"
GBM <- subset(GBM, downsample = 500)


# Rearrange the order of groups.
GBM.list <- SplitObject(GBM, split.by = "infercnv_type_tumor_patient")
GBM.list <- GBM.list[sort(names(GBM.list))]
GBM <- merge(GBM.list[[1]], GBM.list[c(2:25)])


# Get counts matrix.
counts_matrix <- GetAssayData(GBM, slot = "counts")


# Get cell annotations.
cell_Annotations <- GBM@meta.data["infercnv_type_tumor_patient"]
write.table(cell_Annotations, file = paste0(subdir_path, '/cell_Annotations.txt'), quote = F, sep= "\t", col.names = FALSE)


# InferCNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = paste0(subdir_path, '/cell_Annotations.txt'),
                                    gene_order_file = './data/gencode_v19_gen_pos.complete.txt',
                                    delim = "\t",
                                    ref_group_names = c("Macrophage / Microglia", "Oligodendrocyte", "T cell"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = subdir_path, 
                             window_length = win_len[1],
                             cluster_by_groups = TRUE, 
                             denoise = TRUE,
                             HMM = FALSE, 
                             save_rds = FALSE)
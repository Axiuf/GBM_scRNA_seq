# Considering that this is an explorative step to assess the heterogeneity within tumors and different 
# clusters, I provide a test script for convenience of customization and parallelization.


# Set environment
library(Seurat)
library(infercnv)


# Take arguments from the bash script
win_len <- commandArgs(TRUE)
win_len <- as.numeric(win_len)


# Create output dir
data_dir <- "./data/07_GBM_total_merge_InferCNV_test"
dir.create(data_dir)
data_subdir <- paste0(data_dir, "/window_length_", win_len[1])
dir.create(data_subdir)


# Get the Seurat object and sample
GBM <- readRDS("./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
Idents(GBM) <- "infercnv_type_tumor_patient"
GBM <- subset(GBM, downsample = 500)


# Rearrange the order of groups
GBM.list <- SplitObject(GBM, split.by = "infercnv_type_tumor_patient")
GBM.list <- GBM.list[sort(names(GBM.list))]
GBM <- merge(GBM.list[[1]], GBM.list[c(2:25)])


# Get counts matrix
counts_matrix <- GetAssayData(GBM, slot = "counts")


# Get cell annotations
cell_Annotations <- GBM@meta.data["infercnv_type_tumor_patient"]
write.table(cell_Annotations, file = paste0(data_subdir, "/cell_Annotations.txt"),
            quote = F, sep= "\t", col.names = FALSE)


# Infer CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = paste0(data_subdir, "/cell_Annotations.txt"),
                                    gene_order_file = "./data/gencode_v19_gen_pos.complete.txt",
                                    delim = "\t",
                                    ref_group_names = c("Macrophage / Microglia", "Oligodendrocyte", "T cell"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = data_subdir, 
                             window_length = win_len[1],
                             cluster_by_groups = TRUE, 
                             denoise = TRUE,
                             HMM = FALSE, 
                             save_rds = FALSE)
# Set environment
library(Seurat)
library(infercnv)

data_dir <- "./data/08_GBM_total_merge_InferCNV"
dir.create(data_dir)


# Get the Seurat object
GBM <- readRDS("./data/05_GBM_total_merge_filtered_SCT_umap&tsne_Marker.rds")
Idents(GBM) <- "infercnv_type_tumor_patient"


# Rearrange the order of groups
GBM.list <- SplitObject(GBM, split.by = "infercnv_type_tumor_patient")
GBM.list <- GBM.list[sort(names(GBM.list))]
GBM <- merge(GBM.list[[1]], GBM.list[c(2:25)])


# Get counts matrix
counts_matrix <- GetAssayData(GBM, slot = "counts")


# Get cell annotations
cell_Annotations <- GBM@meta.data["infercnv_type_tumor_patient"]
write.table(cell_Annotations, file = paste0(data_dir, "/cell_Annotations.txt"),
            quote = FALSE, sep= "\t", col.names = FALSE)


# Infer CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = paste0(data_dir, "/cell_Annotations.txt"),
                                    gene_order_file = "./data/supplementaries/gencode_v19_gen_pos.complete.txt",
                                    delim = "\t",
                                    ref_group_names = c("Macrophage / Microglia", "Oligodendrocyte", "T cell"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = data_dir, 
                             window_length = 101,
                             cluster_by_groups = TRUE, 
                             denoise = TRUE,
                             HMM = TRUE, 
                             output_format = "pdf",
                             save_rds = TRUE)
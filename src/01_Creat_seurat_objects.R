# Creat seurat objects from the output files of Cell ranger (10 X).
# This script is written specifically for the glioblastoma (GBM) project.
# Make any adjustments necessary for your own data.


library(Seurat)
dir.create("./data/01_Creat_seurat_objects")
dir.create("./data/01_Creat_seurat_objects/excluded")
dir.create("./data/01_Creat_seurat_objects/merged_objects")


# Load data sets and creat Seurat objects individually.
GBM_sample_names <- list.files('./data/Cellranger', pattern = '^P')

GBM_objects <- lapply(GBM_sample_names, function(GBM_sample_name){
  counts_10x <- Read10X(data.dir = paste0("./data/Cellranger/", GBM_sample_name, "/outs/filtered_feature_bc_matrix"))
  raw_object <- CreateSeuratObject(counts = counts_10x, project = GBM_sample_name, min.cells = 3)
  
  saveRDS(raw_object, file = paste0("./data/01_Creat_seurat_objects/", GBM_sample_name, ".rds"))
  
  raw_object
})


# Move P768.rds (in which cell count is too small) to "./data/01_Creat_seurat_objects/excluded".
file.copy("./data/01_Creat_seurat_objects/P768.rds", "./data/01_Creat_seurat_objects/excluded/P768.rds")
file.remove("./data/01_Creat_seurat_objects/P768.rds")

GBM_sample_names <- GBM_sample_names[-13]
GBM_objects <- GBM_objects[-13]


# Merge toal objects.
GBM_total_merge <- merge(GBM_objects[[1]], y = GBM_objects[2:22], add.cell.ids = GBM_sample_names, project = "GBM_total_merge")
saveRDS(GBM_total_merge, file = "./data/01_Creat_seurat_objects/merged_objects/GBM_total_merge.rds")


# Merge objects by patients (in which more than one sample was tested).
P673_merge <- merge(GBM_objects[[7]], y = GBM_objects[[8]], add.cell.ids = GBM_sample_names[7:8], project = "P673_merge")
saveRDS(P673_merge, file = "./data/01_Creat_seurat_objects/merged_objects/P673_merge.rds")

P828_merge <- merge(GBM_objects[[13]], y = GBM_objects[[14]], add.cell.ids = GBM_sample_names[13:14], project = "P828_merge")
saveRDS(P828_merge, file = "./data/01_Creat_seurat_objects/merged_objects/P828_merge.rds")

P843_merge <- merge(GBM_objects[[15]], y = GBM_objects[[16]], add.cell.ids = GBM_sample_names[15:16], project = "P843_merge")
saveRDS(P843_merge, file = "./data/01_Creat_seurat_objects/merged_objects/P843_merge.rds")

P912_merge <- merge(GBM_objects[[20]], y = GBM_objects[[21]], add.cell.ids = GBM_sample_names[20:21], project = "P912_merge")
saveRDS(P912_merge, file = "./data/01_Creat_seurat_objects/merged_objects/P912_merge.rds")

P653_merge <- merge(GBM_objects[[5]], y = GBM_objects[[6]], add.cell.ids = GBM_sample_names[5:6], project = "P653_merge")
saveRDS(P653_merge, file = "./data/01_Creat_seurat_objects/merged_objects/P653_merge.rds")

P679_merge <- merge(GBM_objects[[9]], y = GBM_objects[[10]], add.cell.ids = GBM_sample_names[9:10], project = "P679_merge")
saveRDS(P679_merge, file = "./data/01_Creat_seurat_objects/merged_objects/P679_merge.rds")


# Merge objects by tumor types.
primary_total_merge <- merge(GBM_objects[[1]], y = GBM_objects[c(7, 8, 11, 13:17, 20:22)], add.cell.ids = GBM_sample_names[c(1, 7, 8, 11, 13:17, 20:22)], project = "primary_total_merge")
saveRDS(primary_total_merge, file = "./data/01_Creat_seurat_objects/merged_objects/primary_total_merge.rds")

recurrent_total_merge <- merge(GBM_objects[[2]], y = GBM_objects[c(5, 6, 18, 19)], add.cell.ids = GBM_sample_names[c(2, 5, 6, 18, 19)], project = "recurrent_total_merge")
saveRDS(recurrent_total_merge, file = "./data/01_Creat_seurat_objects/merged_objects/recurrent_total_merge.rds")

Secondary_total_merge <- merge(GBM_objects[[3]], y = GBM_objects[c(4, 9, 10, 12)], add.cell.ids = GBM_sample_names[c(3, 4, 9, 10, 12)], project = "Secondary_total_merge")
saveRDS(Secondary_total_merge, file = "./data/01_Creat_seurat_objects/merged_objects/Secondary_total_merge.rds")

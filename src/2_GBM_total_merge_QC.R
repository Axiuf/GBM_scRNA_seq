# Set environment.
source("requirements.R")
dir.create("./plots/2_GBM_total_merge_QC")
dir.create("./plots/2_GBM_total_merge_QC/percent.mt_VlnPlot")
fig_dpi <- 150

# Get Seurat object.
GBM <- readRDS("./data/raw_objects/merged_objects/GBM_total_merge.rds")


# Add more metrics to meta.data.
GBM[["percent.mt"]] <- PercentageFeatureSet(GBM, pattern = "^MT-")
GBM$log10GenesPerUMI <- log10(GBM$nFeature_RNA) / log10(GBM$nCount_RNA)


# Plot mt_VlnPlot to assess percent.mt by different samples.
for(GBM_sample_name in levels(GBM@active.ident)){
  VlnPlot(GBM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents = GBM_sample_name, pt.size = 0.2, ncol = 3)
  file_name = paste0("percent.mt_VlnPlot_", GBM_sample_name, ".tiff")
  ggsave(filename = file_name, device = "tiff", path = "./plots/2_GBM_total_merge_QC/percent.mt_VlnPlot", dpi = fig_dpi)
}


# Plot percent.mt and nFeature_RNA across nCount_RNA.
plot1 <- FeatureScatter(GBM, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.2, group.by = "orig.ident") + NoLegend()
plot2 <- FeatureScatter(GBM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.2, group.by = "orig.ident") + NoLegend()
plot1 + plot2
ggsave(filename = "percent.mt&nFeature_RNA_by_nCount_FeatureScatter.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, dpi = fig_dpi)


# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/04_SC_quality_control.md
# Plot different metrics to assess the data validity.
metadata <- GBM@meta.data

metadata %>% 
  ggplot(aes(x = orig.ident, fill = orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("nCells")
ggsave(filename = "nCells_by_sample_geom_bar.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 8, dpi = fig_dpi)

metadata %>% 
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell_density") +
  geom_vline(xintercept = 1000)
ggsave(filename = "nCount_RNA_by_sample_geom_density.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, height = 5, dpi = fig_dpi)

metadata %>% 
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell_density") +
  geom_vline(xintercept = 600) 
ggsave(filename = "nFeature_RNA_by_sample_geom_density.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, height = 5, dpi = fig_dpi)

metadata %>% 
  ggplot(aes(x = orig.ident, y=log10(nFeature_RNA), fill = orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells vs nFeature_RNA")
ggsave(filename = "nCells_vs_nFeature_RNA_geom_boxplot.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, height = 5, dpi = fig_dpi)

metadata %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 600) +
  facet_wrap(~orig.ident)
ggsave(filename = "percent.mt_geom_point.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 25, height = 40, dpi = fig_dpi)

metadata %>% 
  ggplot(aes(color = orig.ident, x = percent.mt, fill = orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)
ggsave(filename = "percent.mt_geom_density.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, height = 5, dpi = fig_dpi)

metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave(filename = "log10GenesPerUMI_geom_density.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, height = 5, dpi = fig_dpi)


# Filter cells in the raw data.
GBM_filtered <- subset(x = GBM, subset = (nCount_RNA >= 1000) & 
                            (nFeature_RNA >= 600) & 
                            (log10GenesPerUMI > 0.80) & 
                            (percent.mt < 20))


# Filter genes in the raw data.
counts <- GetAssayData(object = GBM_filtered, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
counts_filtered <- counts[keep_genes, ]
GBM_filtered <- CreateSeuratObject(counts_filtered, meta.data = GBM_filtered@meta.data)
Idents(GBM_filtered) <- "orig.ident"


# plot the filtered data.
plot1 <- FeatureScatter(GBM_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.2, group.by = "orig.ident") + NoLegend()
plot2 <- FeatureScatter(GBM_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.2, group.by = "orig.ident") + NoLegend()
plot1 + plot2
ggsave(filename = "filtered_percent.mt&nFeature_RNA_by_nCount_FeatureScatter.tiff", device = "tiff", path = "./plots/2_GBM_total_merge_QC", width = 14, dpi = fig_dpi)


saveRDS(GBM_filtered, file = "./data/2_GBM_total_merge_filtered.rds")
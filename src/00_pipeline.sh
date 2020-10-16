#!/bin/bash

# Standard work flow
# Rscript ./src/01_Create_seurat_objects.R
# Rscript ./src/02_GBM_total_merge_QC.R
# Rscript ./src/03_1_GBM_total_merge_Normalization1_SCTransform.R
# Rscript ./src/04_GBM_total_merge_Cluster.R
# Rscript ./src/05_GBM_total_merge_Marker.R
# Rscript ./src/08_GBM_total_merge_InferCNV.R # Time consuming
# Rscript ./src/09_2_GBM_tumor_type_SCTIntegrate.R 30
# Rscript ./src/10_2_GBM_tumor_type_Marker.R
# Rscript ./src/11_GBM_tumor_type_Signature_score.R
# Rscript ./src/12_GBM_tumor_type_GO_enrichment.R
# Rscript ./src/13_GBM_non_malignant.R

# Alternative section
# Rscript ./src/03_2_GBM_total_merge_Normalization2_ScaleData.R
# Rscript ./src/06_GBM_total_merge_Simplified_workflow.R
# Rscript ./src/07_GBM_total_merge_InferCNV_test.R 101
# Rscript ./src/09_1_GBM_tumor_type_SCTIntegrate.R 30
# Rscript ./src/09_3_GBM_tumor_type_SCTIntegrate.R 30
# Rscript ./src/10_1_GBM_tumor_type_Marker.R
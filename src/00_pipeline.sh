#!/bin/bash

Rscript ./src/01_Create_seurat_objects.R
Rscript ./src/02_GBM_total_merge_QC.R
Rscript ./src/03_1_GBM_total_merge_Normalization1_SCTransform.R
Rscript ./src/04_GBM_total_merge_Cluster.R
Rscript ./src/05_GBM_total_merge_Marker.R
Rscript ./src/09_2_GBM_tumor_type_SCTIntegrate.R 30
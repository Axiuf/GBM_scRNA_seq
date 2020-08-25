# Import the required packages.
library(dplyr)

library(Seurat)
library(sctransform)

library(patchwork)
library(ggplot2)

library(future)

# Make adjustment by your equipment to run faster.
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 64 * 1024 * 1024^2)
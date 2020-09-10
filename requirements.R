# Import the required packages
library(Seurat)
library(sctransform)

library(dplyr)

library(patchwork)
library(ggplot2)

library(future)
# Some steps in Seurat work flow is time consuming, make adjustment by your equipment to run faster
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 64 * 1024 * 1024^2)

# Change to 300 for officially output
fig_dpi <- 300
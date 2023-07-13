#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)

panc8 <- readRDS(file = "data/panc8_3.0.2.rds")

panc.list <- SplitObject(object = panc8, split.by = "dataset")


panc.list <- lapply(X = panc.list, FUN = SCTransform)

saveRDS(object = panc.list, file = "seurat_objects/panc8.rds")





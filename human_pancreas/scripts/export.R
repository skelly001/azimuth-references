#!/usr/bin/env Rscript
library(Matrix)
library(Seurat) #need seurat >=4.0.0 to work with Azimuth. Need spatstat 1.64-1 for seurat
library(Azimuth)
library(presto)
library(dplyr)
# remotes::install_version("Seurat", version="4.0.0", repos = "http://cran.us.r-project.org", force=T)
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source")
args <- commandArgs(trailingOnly = TRUE)

ref.dir <- "reference/"
ob.dir <- "seurat_objects/"
ref <- readRDS(file = "seurat_objects/pancreas_integrated.rds")
annotations <- readRDS(file = "seurat_objects/annotations.rds")
Idents(object = ref) <- annotations

if ("remove" %in% levels(x = ref)) {
  ref <- subset(x = ref, idents = "remove", invert = TRUE)
  ref <- RunPCA(object = ref, verbose = FALSE)
}
ref$annotation.l1 <- Idents(object = ref)
ref <- RunUMAP(object = ref, dims = 1:30, return.model = TRUE)
full.ref <- ref

colormap <- list(annotation.l1 = CreateColorMap(object = ref, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

ref2 <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "integrated",
  metadata = c("annotation.l1"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))
saveRDS(object = full.ref, file = file.path(ob.dir, "fullref.Rds"))


library(Matrix)
library(Seurat) #need seurat >=4.0.0 to work with Azimuth. Need spatstat 1.64-1 for seurat
library(Azimuth)
library(presto)
library(dplyr)


final_ref <- readRDS(file = "seurat_objects/final_ref.rds")

# Meta data
meta <- final_ref@meta.data

# Count data
sct_counts <- as.matrix(final_ref@assays[["SCT"]]@counts[,1:50]) %>% 
  subset.matrix(rownames(.) %in% c("INS", "GCG"))

sct_data <- as.matrix(final_ref@assays[["SCT"]]@data[,1:50])%>% 
  subset.matrix(rownames(.) %in% c("INS", "GCG"))

identical(sct_counts, sct_data)

# Add meta to counts








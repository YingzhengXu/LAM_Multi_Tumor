library(Seurat)
library(SeuratDisk)
library(anndata)

# Read all H5AD files in a folder
list.h5ad <- list()
h5_files <- list.files("h5ad", pattern = "\\.h5ad$")
for (i in seq_along(h5_files)) {
  file_path <- file.path("h5ad", h5_files[i])
  list.h5ad[[i]] <- anndata::read_h5ad(file_path)
}
names(list.h5ad) <- h5_files

# Convert to Seurat objects
list.seu <- list()
for (i in seq_along(list.h5ad)) {
  cat("Creating Seurat object for", names(list.h5ad)[i], "\n")
  seu_obj <- CreateSeuratObject(counts = list.h5ad[[i]]$transpose()$X,
                                meta.data = list.h5ad[[i]]$obs,
                                project = names(list.h5ad)[i])
  list.seu[[i]] <- seu_obj
}
names(list.seu) <- h5_files

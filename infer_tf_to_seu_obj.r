library(Seurat)

# Directory where AUCell files are saved
aucell_dir <- "output/"

# List AUCell files (mac example, adapt as needed)
aucell_files <- list.files(aucell_dir, pattern = "^mac_AUCell_.*\\.csv$", full.names = TRUE)

# Load your Seurat objects (assume stored in RDS)
seurat_files <- list.files("seurat_rds/", pattern = "\\.rds$", full.names = TRUE)
names(seurat_files) <- sub("\\.rds$", "", basename(seurat_files))

# Loop through AUCell files
seurat_list <- list()
for (file in aucell_files) {
  # Extract tumor name from filename
  tumor_name <- sub("^mac_AUCell_(.*)\\.csv$", "\\1", basename(file))
  
  # Load AUCell
  aucell_scores <- read.csv(file, row.names = 1)
  aucell_scores <- as.matrix(aucell_scores)
  
  # Load Seurat
  seu_path <- seurat_files[[tumor_name]]
  if (is.null(seu_path)) {
    message(paste("No Seurat object for", tumor_name))
    next
  }
  seu <- readRDS(seu_path)
  
  # Align cells
  common_cells <- intersect(colnames(seu), colnames(aucell_scores))
  if (length(common_cells) == 0) {
    message(paste("No overlapping cells for", tumor_name))
    next
  }
  aucell_aligned <- aucell_scores[, common_cells]
  
  # Add AUCell as assay
  seu[["TF"]] <- CreateAssayObject(counts = aucell_aligned)
  DefaultAssay(seu) <- "TF"
  
  # Save updated Seurat object
  saveRDS(seu, file = paste0("seurat_rds_with_TF/", tumor_name, "_TF.rds"))
  
  # Optionally store for plotting
  seurat_list[[tumor_name]] <- seu
}

# Example visualization for one or all
FeaturePlot(seurat_list[["Breast_BC"]], features = c("TREM2", "SPP1"))
DotPlot(seurat_list[["Breast_BC"]], features = c("TREM2", "SPP1")) + coord_flip()

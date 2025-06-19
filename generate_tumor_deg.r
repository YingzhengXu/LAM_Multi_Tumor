library(Seurat)
library(anndata)

# 1. Read all h5ad files & convert to seu obj
h5ad_files <- list.files("h5ad", full.names = TRUE)
list.seu <- list()

for (file in h5ad_files) {
  h5 <- anndata::read_h5ad(file)
  seu <- CreateSeuratObject(counts = h5$transpose()$X, meta.data = h5$obs, project = basename(file))
  seu$dataID <- basename(file)
  list.seu[[basename(file)]] <- seu
}


# 2. Subset to macrophages / monocytes / LAM relevant cells
for (i in seq_along(list.seu)) {
  list.seu[[i]] <- subset(list.seu[[i]], annotation_immune %in% c("Macrophage", "Monocyte"))
}


# 3. DE analysis Tumor vs control
for (i in seq_along(list.seu)) {
  Idents(list.seu[[i]]) <- list.seu[[i]]$condition
}

macdeg <- list()
for (i in seq_along(list.seu)) {
  macdeg[[names(list.seu)[i]]] <- FindMarkers(list.seu[[i]], ident.1 = "Tumor", ident.2 = "Normal", max.cells.per.ident = 500)
}

# 4. Generate DEG expression matrices
dir.create("macdegmtx", showWarnings = FALSE)
for (i in seq_along(list.seu)) {
  genes <- rownames(macdeg[[i]])
  expr_matrix <- GetAssayData(list.seu[[i]], slot = "data")[genes, ]
  expr_df <- as.data.frame(as.matrix(expr_matrix))
  write.csv(expr_df, file = paste0("macdegmtx/", names(list.seu)[i], ".csv"), row.names = TRUE)
}


# 5. Generate LAM-specific DEG matrices
lamdeg <- list()
for (i in seq_along(list.seu)) {
  lam_cells <- subset(list.seu[[i]], CLUSTER == "TREM2+")  # Example: LAM cluster
  other_cells <- subset(list.seu[[i]], CLUSTER != "TREM2+")
  
  lamdeg[[names(list.seu)[i]]] <- FindMarkers(lam_cells, ident.1 = "LAM", ident.2 = "Other", max.cells.per.ident = 500)
}

dir.create("lamdegmtx", showWarnings = FALSE)
for (i in seq_along(list.seu)) {
  genes <- rownames(lamdeg[[i]])
  expr_matrix <- GetAssayData(list.seu[[i]], slot = "data")[genes, ]
  expr_df <- as.data.frame(as.matrix(expr_matrix))
  write.csv(expr_df, file = paste0("lamdegmtx/", names(list.seu)[i], ".csv"), row.names = TRUE)
}



# preprocess, integrate, cluster
for (i in seq_along(list.seu)) {
  list.seu[[i]] <- list.seu[[i]] %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.4)
}

# find DEGs between two conditions within a dataset
deg_results <- FindMarkers(list.seu[[1]], ident.1 = "Tumor", ident.2 = "Normal")
write.csv(deg_results, paste0("deg_", names(list.seu)[1], ".csv"))

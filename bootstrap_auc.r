# Load required libraries
library(dplyr)

# Directory containing AUCell scores
aucell_dir <- "output"  # adjust to your actual output path
aucell_files <- list.files(aucell_dir, pattern = "mac_AUCell_.*\\.csv$", full.names = TRUE)

# Load all AUCell matrices
aucell_matrices <- lapply(aucell_files, function(f) {
  mat <- read.csv(f, row.names = 1)
  mat
})
names(aucell_matrices) <- gsub("mac_AUCell_|\\.csv", "", basename(aucell_files))

# Bootstrap function
bootstrap_auc <- function(auc_matrix, n_bootstrap = 1000) {
  tf_ranks <- matrix(NA, nrow = n_bootstrap, ncol = nrow(auc_matrix))
  rownames(tf_ranks) <- 1:n_bootstrap
  colnames(tf_ranks) <- rownames(auc_matrix)
  
  set.seed(123)
  for (i in 1:n_bootstrap) {
    sampled_cols <- sample(ncol(auc_matrix), replace = TRUE)
    resampled <- auc_matrix[, sampled_cols, drop = FALSE]
    tf_scores <- rowSums(resampled)
    tf_ranks[i, ] <- rank(-tf_scores)  # rank: smaller = stronger
  }
  avg_ranks <- colMeans(tf_ranks)
  stable_rank <- sort(avg_ranks)
  return(stable_rank)
}

# Apply bootstrap to all AUCell matrices
bootstrap_results <- lapply(aucell_matrices, bootstrap_auc)

# Example: view top 10 TFs for one tumor type
print(bootstrap_results[["Breast_BC"]][1:10])

# Optionally save bootstrap results
for (name in names(bootstrap_results)) {
  write.csv(bootstrap_results[[name]], file = paste0("output/bootstrap_rank_", name, ".csv"))
}

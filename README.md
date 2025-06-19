# Pan-Cancer Macrophage GRN and Survival Analysis

This repository contains code for a complete pipeline to:
- Build Seurat objects from `.h5ad` files
- Identify macrophage and LAM DEGs across tumor types
- Run GRNBoost2 and calculate AUCell TF activity scores
- Integrate AUCell scores into Seurat and visualize TF activity
- Perform bootstrap analysis of TF ranks
- Conduct survival analysis of LAM features using TCGA data

---

## Pipeline Overview

### 1. **Generate DEGs**
`generate_tumor_deg.r`  
- Reads `.h5ad` files
- Creates Seurat objects
- Identifies DEGs for macrophages (macdegmtx/) and LAM (lamdegmtx/)
- Exports DEG expression matrices

---

### 2. **GRN + AUCell**
`grn_aucell.py`  
- Loads DEG matrices
- Infers GRNs using GRNBoost2
- Saves networks (`mac_network/`, `lam_network/`)
- Computes AUCell TF scores
- Saves AUCell scores (`output/`)

---

### 3. **Integrate AUCell scores**
`infer_tf_to_seu_obj.r`  
- Reads AUCell scores
- Adds as new assay to Seurat objects
- Visualizes with `FeaturePlot`, `DotPlot`

---

### 4. **Bootstrap TF rankings**
`bootstrap_auc.r`  
- Performs resampling on AUCell matrices
- Ranks TFs based on stability across bootstraps

---

### 5. **Survival analysis**
`tcga_survival_of_lam_features.r`  
- Analyzes association between LAM signature / TF activity and patient survival using TCGA data

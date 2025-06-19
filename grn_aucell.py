import os
import pandas as pd
from arboreto.algo import grnboost2
from pyscenic.aucell import aucell
from pyscenic.rss import Regulon
from dask.diagnostics import ProgressBar

# === Directories ===
BASE_DIR = "/Users/yingzhengxu/Desktop/pan cancer"
MAC_DEG_DIR = os.path.join(BASE_DIR, "macdegmtx")
LAM_DEG_DIR = os.path.join(BASE_DIR, "lamdegmtx")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
TF_FILE = os.path.join(BASE_DIR, "tf.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# === Load transcription factors ===
tf_names = pd.read_csv(TF_FILE, index_col=0, header=None).index.tolist()

# === Helper functions ===
def load_expression_matrices(directory, transpose=True):
    """Load CSV matrices and optionally transpose"""
    matrices = []
    filenames = sorted([f for f in os.listdir(directory) if f.endswith(".csv")])
    for f in filenames:
        df = pd.read_csv(os.path.join(directory, f), index_col=0)
        if transpose:
            df = df.T
        matrices.append(df)
    return matrices, filenames

def infer_grnboost2(data_list, tf_names):
    """Run GRNBoost2 on a list of expression matrices"""
    networks = []
    with ProgressBar():
        for i, data in enumerate(data_list):
            print(f"[GRNBoost2] Inferring network {i + 1}/{len(data_list)}")
            net = grnboost2(expression_data=data, tf_names=tf_names)
            networks.append(net)
    return networks

def save_networks(networks, filenames, prefix):
    """Save networks to CSV files"""
    for net, fname in zip(networks, filenames):
        out_path = os.path.join(OUTPUT_DIR, f"{prefix}_network_{fname}")
        net.to_csv(out_path, index=False)
        print(f"[Saved] {out_path}")

def convert_to_regulons(network_df):
    """Convert GRNBoost2 output to regulons for AUCell"""
    regulons = []
    for tf, group in network_df.groupby("TF"):
        gene2weight = dict(zip(group["target"], group["importance"]))
        regulons.append(Regulon(name=tf, gene2weight=gene2weight))
    return regulons

def calculate_aucell_scores(expr_df, network_df):
    """Calculate AUCell scores from expression data + network"""
    regulons = convert_to_regulons(network_df)
    return aucell(expr_df, regulons)

def save_aucell_scores(data_list, networks, filenames, prefix):
    """Compute AUCell scores and save them"""
    for df, net, fname in zip(data_list, networks, filenames):
        print(f"[AUCell] Computing for {fname}")
        regulons = convert_to_regulons(net)
        scores = aucell(df, regulons)
        out_path = os.path.join(OUTPUT_DIR, f"{prefix}_AUCell_{fname}")
        scores.to_csv(out_path)
        print(f"[Saved AUCell] {out_path}")

# === Pipeline ===
# Load data
mac_data, mac_files = load_expression_matrices(MAC_DEG_DIR)
lam_data, lam_files = load_expression_matrices(LAM_DEG_DIR)

# GRNBoost2
mac_networks = infer_grnboost2(mac_data, tf_names)
lam_networks = infer_grnboost2(lam_data, tf_names)

save_networks(mac_networks, mac_files, "mac")
save_networks(lam_networks, lam_files, "lam")

# Run AUCell and save
save_aucell_scores(mac_data, mac_networks, mac_files, "mac")
save_aucell_scores(lam_data, lam_networks, lam_files, "lam")


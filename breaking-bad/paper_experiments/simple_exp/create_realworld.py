import argparse
import pandas as pd
import numpy as np
import os
import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

robjects.r('options(warn=-1)')

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", required=True,
                    help="Output dataset name (e.g. ARTH150, ECOLI70, MAGIC-IRRI, MAGIC-NIAB).")
parser.add_argument("--rds", required=True,
                    help="Filename of the .rds file located next to this script (e.g. arth150.rds).")
parser.add_argument("--samples", type=int, default=10000,
                    help="Number of samples per generated dataset (default: 10000).")
parser.add_argument("--runs", type=int, default=30,
                    help="Number of seeds, generates data for seeds 1..runs (default: 30).")
args = parser.parse_args()

DATASET_NAME = args.dataset
N_SAMPLES = args.samples
SEEDS = list(range(1, args.runs + 1))

user_lib = os.path.expanduser("~/R_libs")
if not os.path.exists(user_lib):
    os.makedirs(user_lib)

robjects.r(f'.libPaths(c("{user_lib}", .libPaths()))')

utils = importr('utils')

try:
    bnlearn = importr('bnlearn')
    print("bnlearn is loaded successfully.")
except Exception:
    print(f"bnlearn not found. Installing into {user_lib}...")
    package_url = "https://cran.r-project.org/src/contrib/Archive/bnlearn/bnlearn_4.8.1.tar.gz"
    r_null = robjects.r('NULL')
    utils.install_packages(package_url, repos=r_null, type="source", lib=user_lib)
    bnlearn = importr('bnlearn')
    print("bnlearn installed and loaded successfully.")

current_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(current_dir, "..", "..", "data", "Realworld", DATASET_NAME))
os.makedirs(output_dir, exist_ok=True)

rds_path = os.path.join(current_dir, args.rds)

if not os.path.exists(rds_path):
    print(f"\n[Error] '{rds_path}' files does not exist.")
    sys.exit(1)

print(f"Loading local network file: {rds_path}")
r_bn_object = robjects.r['readRDS'](rds_path)

print(f"\n[{DATASET_NAME}] Extracting Adjacency Matrix...")

r_amat = bnlearn.amat(r_bn_object)

with localconverter(robjects.default_converter + pandas2ri.converter):
    adj_data = robjects.conversion.rpy2py(r_amat)

r_nodes = bnlearn.nodes(r_bn_object)
node_names = list(r_nodes)
adj_df = pd.DataFrame(np.array(adj_data), index=node_names, columns=node_names)

adj_csv_name = f"True_DAG_{DATASET_NAME}.csv"
adj_csv_path = os.path.join(output_dir, adj_csv_name)
adj_df.to_csv(adj_csv_path, index=False) 

adj_npy_name = f"True_DAG_{DATASET_NAME}.npy"
adj_npy_path = os.path.join(output_dir, adj_npy_name)
adj_numpy_c_order = np.ascontiguousarray(adj_df.values)
np.save(adj_npy_path, adj_numpy_c_order)

print(f" - Saved: {adj_csv_name}")
print(f" - Saved: {adj_npy_name}")

print(f"\nGenerating {len(SEEDS)} datasets with {N_SAMPLES} samples each...")

for seed in SEEDS:
    robjects.r['set.seed'](seed)
    
    r_data = bnlearn.rbn(r_bn_object, n=N_SAMPLES)
    
    with localconverter(robjects.default_converter + pandas2ri.converter):
        df = robjects.conversion.rpy2py(r_data)
    
    file_base = f"Data_{DATASET_NAME}_{N_SAMPLES}_{seed}"
    
    csv_name = f"{file_base}.csv"
    csv_path = os.path.join(output_dir, csv_name)
    df.to_csv(csv_path, index=False)
    
    npy_name = f"{file_base}.npy"
    npy_path = os.path.join(output_dir, npy_name)
    data_numpy_c_order = np.ascontiguousarray(df.values)
    np.save(npy_path, data_numpy_c_order)
    
    print(f"[{seed}/{len(SEEDS)}] Saved: {csv_name}")

print(f"\nAll Done! Files are in '{output_dir}' folder.")
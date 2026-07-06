import os
import argparse
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

robjects.r('options(warn=-1)')

DATASET_NAME = "Sachs"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="sachs.data.txt")
    args = parser.parse_args()

    user_lib = os.path.expanduser("~/R_libs")
    os.makedirs(user_lib, exist_ok=True)
    robjects.r(f'.libPaths(c("{user_lib}", .libPaths()))')
    utils = importr('utils')
    try:
        bnlearn = importr('bnlearn')
    except Exception:
        package_url = "https://cran.r-project.org/src/contrib/Archive/bnlearn/bnlearn_4.8.1.tar.gz"
        utils.install_packages(package_url, repos=robjects.r('NULL'), type="source", lib=user_lib)
        bnlearn = importr('bnlearn')

    current_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.abspath(os.path.join(current_dir, "..", "..", "data", "Realworld", DATASET_NAME))
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output Directory: {output_dir}")
    data_path = os.path.join(current_dir, args.data)
    if not os.path.exists(data_path):
        print(f"❌ File not found: {data_path}")
        return

    df = pd.read_csv(data_path, sep="\t")
    node_names = list(df.columns)
    N = len(df)
    print(f"Loaded {N} rows × {len(node_names)} columns: {node_names}")

    file_base = f"Data_{DATASET_NAME}_{N}_1"
    df.to_csv(os.path.join(output_dir, f"{file_base}.csv"), index=False)
    np.save(os.path.join(output_dir, f"{file_base}.npy"), np.ascontiguousarray(df.values))
    print(f"✅ Saved Data: {file_base}.csv / .npy")

    sachs_dag_str = (
        "[PKC][PKA|PKC][Raf|PKC:PKA][Mek|PKC:PKA:Raf]"
        "[Erk|Mek:PKA][Akt|Erk:PKA][P38|PKC:PKA][Jnk|PKC:PKA]"
        "[Plcg][PIP3|Plcg][PIP2|Plcg:PIP3]"
    )
    r_dag = robjects.r['model2network'](sachs_dag_str)
    r_amat = bnlearn.amat(r_dag)
    
    with localconverter(robjects.default_converter + pandas2ri.converter):
        adj_data = robjects.conversion.rpy2py(r_amat)
        
    dag_nodes = list(bnlearn.nodes(r_dag))
    adj_df = pd.DataFrame(np.array(adj_data), index=dag_nodes, columns=dag_nodes)

    adj_df = adj_df.reindex(index=node_names, columns=node_names)
    
    gt_csv_name = f"True_DAG_{DATASET_NAME}.csv"
    adj_df.to_csv(os.path.join(output_dir, gt_csv_name), index=False)
    np.save(os.path.join(output_dir, f"True_DAG_{DATASET_NAME}.npy"), np.ascontiguousarray(adj_df.values))
    print(f"✅ Saved True DAG: {gt_csv_name}")

if __name__ == "__main__":
    main()
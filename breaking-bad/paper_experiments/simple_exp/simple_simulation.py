import numpy as np
import sys
import os
sys.path.append("..")

from utils_lges import dag_to_cpdag, get_random_data, RANDOM_SEED
from sempler.generators import dag_avg_deg
import networkx as nx
import matplotlib.pyplot as plt


# =========================
# Data generation function
# =========================
def simulation_dag_data(
    n_samples,
    p,
    avg_deg,
    given_graph=None,
    graph_seed=RANDOM_SEED,
    data_seed=None
):
    # ---- generate / load DAG (structure seed fixed) ----
    if given_graph is None:
        G_true = dag_avg_deg(
            p, avg_deg, 1, 1, random_state=graph_seed
        )
    else:
        G_true = np.array(given_graph)

    # ---- generate data (data seed depends on sample size) ----
    if data_seed is None:
        data_seed = graph_seed + n_samples

    data = get_random_data(
        G_true,
        n_samples,
        model="linear-gaussian",
        random_state=data_seed
    )[0]

    return G_true, data


# =========================
# Main
# =========================
if __name__ == "__main__":

    # -------------------------
    # Experiment grid
    # -------------------------
    # You can adjust these lists as needed
    # ...existing code...
    # Experiment grid
    # -------------------------
    # You can adjust these lists as needed
    p_list = [20] 
    n_samples_list = [1000]
    avg_deg_list = [2.0]
    
    # Number of graph seeds (1 to 10)
    n_graph_seeds = 30

    # -------------------------
    # Paths
    # -------------------------
    BASE_DIR = "../../data"
    GT_DIR = os.path.join(BASE_DIR, "ground_truth")
    DATA_DIR = os.path.join(BASE_DIR, "samples")

    os.makedirs(GT_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)

    print(f"Generating data for p={p_list}, deg={avg_deg_list}, samples={n_samples_list}, seeds=1..{n_graph_seeds}")

    # -------------------------
    # Generate datasets
    # -------------------------
    for p in p_list:
        for avg_deg in avg_deg_list:
            for g_num in range(1, n_graph_seeds + 1):
                
                # 1. DAG generation (graph_seed = g_num ensures different structure per g_num)
                G_true, _ = simulation_dag_data(
                    n_samples=1,          # dummy
                    p=p,
                    avg_deg=avg_deg,
                    given_graph=None,
                    graph_seed=g_num      # Seed 1~10
                )

                # --- Save Ground Truth (NPY) ---
                gt_filename_npy = f"True_DAG_var{p}_avg_deg{avg_deg}_graph_num{g_num}.npy"
                gt_path_npy = os.path.join(GT_DIR, gt_filename_npy)
                np.save(gt_path_npy, G_true)

                # --- Save Ground Truth (CSV with header) ---
                gt_filename_csv = f"True_DAG_var{p}_avg_deg{avg_deg}_graph_num{g_num}.csv"
                gt_path_csv = os.path.join(GT_DIR, gt_filename_csv)
                
                # Create header string: X0,X1,...,Xp-1
                header_str = ",".join([f"X{i}" for i in range(p)])
                # Save as integer (fmt="%d") with header, removing the hash mark (comments="")
                np.savetxt(gt_path_csv, G_true, delimiter=",", fmt="%d", header=header_str, comments="")

                print(f"[Graph {g_num}] Saved True DAG -> {gt_filename_csv} & .npy (deg={avg_deg})")

                # 2. Generate Data for each sample size
                for n_samples in n_samples_list:
                    # Data seed ensures reproducibility but varies with sample size
                    current_data_seed = g_num * 100000 + n_samples
                    
                    _, data = simulation_dag_data(
                        n_samples=n_samples,
                        p=p,
                        avg_deg=avg_deg,
                        given_graph=G_true,   # Fix the DAG structure
                        data_seed=current_data_seed
                    )

                    # Save Data (NPY)
                    data_filename = (
                        f"Data_var{p}_avg_deg{avg_deg}_n_samples{n_samples}_graph_num{g_num}.npy"
                    )
                    data_path = os.path.join(DATA_DIR, data_filename)
                    np.save(data_path, data)

                    # print(f"    -> Saved Data (n={n_samples})")

    print("✅ All datasets generated successfully.")
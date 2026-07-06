import numpy as np
import sys
import os
sys.path.append("..")

from utils_lges import dag_to_cpdag, get_random_data, RANDOM_SEED
from sempler.generators import dag_avg_deg

def simulation_dag_data(n_samples, p, avg_deg, model_type="non-linear-anm", given_graph=None, graph_seed=RANDOM_SEED, data_seed=None):
    if given_graph is None:
        G_true = dag_avg_deg(p, avg_deg, 1, 1, random_state=graph_seed)
    else:
        G_true = np.array(given_graph)

    if data_seed is None:
        data_seed = graph_seed + n_samples

    data = get_random_data(
        G_true,
        n_samples,
        model=model_type,
        random_state=data_seed
    )
    
    # get_random_data
    if isinstance(data, list):
        data = data[0]

    return G_true, data

if __name__ == "__main__":
    p_list = [20] 
    n_samples_list = [10000]
    avg_deg_list = [5.0]
    n_graph_seeds = 30
    
    models = [("non-linear-anm", "nonlinear"), ("linear-gumbel", "gumbel")]

    BASE_DIR = "../../data"
    GT_DIR = os.path.join(BASE_DIR, "ground_truth")
    DATA_DIR = os.path.join(BASE_DIR, "samples")

    os.makedirs(GT_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)

    print(f"Generating data for p={p_list}, deg={avg_deg_list}, samples={n_samples_list}, seeds=1..{n_graph_seeds}")

    for model_name, prefix in models:
        print(f"\n🚀 Starting {model_name} generation...")
        for p in p_list:
            for avg_deg in avg_deg_list:
                for g_num in range(1, n_graph_seeds + 1):
                    
                    G_true, _ = simulation_dag_data(
                        n_samples=1,
                        p=p,
                        avg_deg=avg_deg,
                        model_type=model_name,
                        given_graph=None,
                        graph_seed=g_num
                    )

                    # True DAG
                    gt_filename_npy = f"True_DAG_{prefix}_var{p}_avg_deg{avg_deg}_graph_num{g_num}.npy"
                    gt_path_npy = os.path.join(GT_DIR, gt_filename_npy)
                    np.save(gt_path_npy, G_true)

                    gt_filename_csv = f"True_DAG_{prefix}_var{p}_avg_deg{avg_deg}_graph_num{g_num}.csv"
                    gt_path_csv = os.path.join(GT_DIR, gt_filename_csv)
                    header_str = ",".join([f"X{i}" for i in range(p)])
                    np.savetxt(gt_path_csv, G_true, delimiter=",", fmt="%d", header=header_str, comments="")

                    for n_samples in n_samples_list:
                        current_data_seed = g_num * 100000 + n_samples
                        
                        _, data = simulation_dag_data(
                            n_samples=n_samples,
                            p=p,
                            avg_deg=avg_deg,
                            model_type=model_name,
                            given_graph=G_true,
                            data_seed=current_data_seed
                        )

                        # Data
                        data_filename = f"Data_{prefix}_var{p}_avg_deg{avg_deg}_n_samples{n_samples}_graph_num{g_num}.npy"
                        data_path = os.path.join(DATA_DIR, data_filename)
                        np.save(data_path, data)

                        print(f"[Graph {g_num} | {prefix}] Saved Data -> {data_filename}")

    print("✅ All datasets generated successfully.")
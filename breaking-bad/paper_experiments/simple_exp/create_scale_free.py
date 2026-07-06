import networkx as nx
import numpy as np
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__)) # simple_exp
parent_dir = os.path.dirname(current_dir) # paper_experiments

if parent_dir not in sys.path:
    sys.path.append(parent_dir)

try:
    from utils_lges import get_random_data, RANDOM_SEED
except ImportError:
    print("Warning: Could not import utils_lges from parent directory. Trying '../'")
    sys.path.append("..")
    from utils_lges import get_random_data, RANDOM_SEED

def create_scale_free_DAG(n, alpha=0.41, beta=0.54, gamma=0.05, delta_in=0.2, delta_out=0, seed=42):

    sc_graph = nx.scale_free_graph(
        n=n, 
        alpha=alpha, 
        beta=beta, 
        gamma=gamma, 
        delta_in=delta_in, 
        delta_out=delta_out, 
        seed=seed
    )

    # (MultiDiGraph -> DiGraph)
    G = nx.DiGraph(sc_graph)

    G.remove_edges_from(list(nx.selfloop_edges(G)))
    
    DAG_edges = [(u, v) for u, v in G.edges() if u > v]

    DAG = nx.DiGraph()
    DAG.add_nodes_from(range(n)) 
    DAG.add_edges_from(DAG_edges)

    return DAG

def create_scale_free_DAG_refined(n,m=4,seed=42):
    G_undirected=nx.barabasi_albert_graph(n,m,seed)
    nodes=list(G_undirected.nodes())
    np.random.shuffle(nodes)
    order_map={node: i for i,node in enumerate(nodes)}

    G=nx.DiGraph()
    G.add_nodes_from(G_undirected.nodes())

    for u, v in G_undirected.edges():
        if order_map[u]<order_map[v]:  
            G.add_edge(u, v) 
        else:
            G.add_edge(v, u) 
    
    threshold=3
    hub_nodes = [n for n, d in G.in_degree() if d >= threshold]
    count = len(hub_nodes)
    print(f"Hub Node: {count}")
    
    return G

# =========================
# Main Execution
# =========================
if __name__ == "__main__":

    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    BASE_DIR = os.path.join(os.path.dirname(base_path), "data")
    GT_DIR = os.path.join(BASE_DIR, "ground_truth")
    DATA_DIR = os.path.join(BASE_DIR, "samples")

    os.makedirs(GT_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)

    var_list = [50]
    samples_list = [1000]
    seeds = 30

    m_value = 4

    print(f"Creating Scale-Free Datasets (m={m_value})...")

    for num_var in var_list:
        for seed in range(1, seeds+1):

            my_dag = create_scale_free_DAG_refined(
                n=num_var,
                m=m_value,
                seed=seed
            )

            dag_adjacency = nx.to_numpy_array(
                my_dag,
                nodelist=sorted(my_dag.nodes()),
                weight=None
            )

            # ---------- Ground Truth ----------
            gt_filename_csv = f"ScaleFree_m{m_value}_DAG_var{num_var}_seed{seed}.csv"
            gt_path_csv = os.path.join(GT_DIR, gt_filename_csv)

            header_str = ",".join([f"X{i}" for i in range(num_var)])
            np.savetxt(gt_path_csv, dag_adjacency, delimiter=",",
                       fmt="%d", header=header_str, comments="")

            np.save(gt_path_csv.replace(".csv", ".npy"), dag_adjacency)

            print(f"[Var{num_var}-Seed{seed}] DAG Saved (m={m_value}). Edges: {my_dag.number_of_edges()}")

            # ---------- Data ----------
            for num_samples in samples_list:

                generated_result = get_random_data(
                    dag_adjacency,
                    num_samples,
                    model="linear-gaussian",
                    random_state=seed * 10000 + num_samples
                )

                data = generated_result[0] if isinstance(generated_result, (list, tuple)) else generated_result

                data_filename = f"ScaleFree_m{m_value}_Data_var{num_var}_n{num_samples}_seed{seed}.npy"
                data_path = os.path.join(DATA_DIR, data_filename)

                np.save(data_path, data)

    print("✅ All Scale-Free datasets generated.")
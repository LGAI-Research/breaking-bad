#!/usr/bin/env python3
import argparse
import json
import os
import time
from datetime import datetime
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from pdag import PDAG, compute_shd_cpdag, pdag_to_dag
from utils_lges import cpdag_shd, dag_to_cpdag
import warnings

warnings.filterwarnings("ignore")



def order_divergence(order, adj):
    err = 0
    for i in range(len(order)):
        err += adj[order[i+1:], order[i]].sum()
    return err



def parse_args():
    p = argparse.ArgumentParser(
        description="Load search result.csv, compute SHD/ecSHD, and save results.json.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )


    # result.csv
    p.add_argument("--result", type=str, required=True,
                   help="Path to search result adjacency CSV (CPDAG encoded)")
    p.add_argument("--stats", type=str, required=True,
                   help="Path to search stats CSV")
    
    #ground_truth
    p.add_argument("--ground_truth", type=str, required=True,
                   help="Path to ground truth adjacency CSV (DAG encoded)")

    default_tag = datetime.now().strftime("run-%Y%m%d-%H%M%S")
    p.add_argument("--tag", type=str, default=default_tag,
                   help="Run tag used to create output subdirectory")
    p.add_argument("--out-root", type=str, default="./output",
                   help="Root directory to store outputs (will create <out-root>/<tag>)")

    return p.parse_args()


def load_square_csv(path: Path) -> np.ndarray:
    try:
        M = pd.read_csv(path, header=None).values
        if M.ndim == 2 and M.shape[0] == M.shape[1]:
            return M
    except Exception:
        pass
    # fallback: header=0
    M = pd.read_csv(path, header=0).values
    if M.ndim == 2 and M.shape[0] == M.shape[1]:
        return M
    raise ValueError(f"result not square: {path} got {M.shape}")


def main():
    args = parse_args()
    gt_path = os.path.join(args.ground_truth)
    result_path = Path(args.result).resolve()

    if not os.path.isfile(gt_path):
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")
    if not result_path.exists():
        raise FileNotFoundError(f"Result CSV not found: {result_path}")

    # ----- (CPDAG) -----
    G_pred = nx.DiGraph(pd.read_csv(result_path, header=0).values)
    G_pred_np = np.array(pd.read_csv(result_path, header=0).values)
    stats = pd.read_csv(args.stats, index_col=0, header=None).T.to_dict("records")[0]
    res = dict()
    for k, v in stats.items():
        res[k] = v


    # ----- (DAG) -----
    G_true = nx.DiGraph(pd.read_csv(gt_path, header=0).values)
    G_true_np = np.array(pd.read_csv(gt_path, header=0).values)
    if not nx.is_directed_acyclic_graph(G_true):
        raise ValueError("G_true is not a DAG (cycle detected). Check your ground-truth file.")

    #check Graph is a DAG
    print("G_true is a DAG: ", nx.is_directed_acyclic_graph(G_true)) # True
    print("G_pred is a DAG: ", nx.is_directed_acyclic_graph(G_pred)) # False -> CPDAG

    try:
        shd = PDAG.from_digraph(G_pred).shd_against_dag(G_true) 
    except Exception as e:
        print(f"[Warn] Standard SHD calc failed ({e}).")
        shd = -1

    # 2. ecSHD (Equivalence Class SHD)
    # Note: It is reported in the paper main result as SHD !

    cpdag_true_lges = dag_to_cpdag(G_true_np)
    try:
        ecshd_res = cpdag_shd(cpdag_true_lges, G_pred_np)
        ecshd = ecshd_res[0] if isinstance(ecshd_res, tuple) else ecshd_res
    except Exception as e:
        print(f"[Warn] Numpy ecSHD calc failed ({e}).")
        ecshd = -1

    res["shd"] = shd
    res["ecshd"] = ecshd

    # 3. Order Divergence
    try:
        pi_pred = list(nx.topological_sort(pdag_to_dag(PDAG.from_digraph(G_pred))))
        OD = order_divergence(pi_pred, G_true_np)
    except ValueError:
        # "No valid node x found..." or "Cycle detected"
        OD = -1
        print("[Warn] Could not compute Order Divergence (Graph is inconsistent/cyclic).")
    
    res["order_divergence"] = OD


    print(f"Loaded result: {result_path.name}")
    print(f"SHD (DAG vs PDAG): {shd}")
    print(f"ecSHD (Direct Comparison): {ecshd}")
    print(f"Order divergence: {OD}")


    # -- sanity check using utils from LGES--
    cpdag_true_lges = dag_to_cpdag(G_true_np)
    ecshd_lges = cpdag_shd(cpdag_true_lges, G_pred_np)
    print("ecSHD (equivalence-class) from LGES git: ", ecshd_lges)
    print("Number of edges in G_true: ", PDAG.from_digraph(G_true).number_of_edges)
    print("Number of edges in G_pred: ", PDAG.from_digraph(G_pred).number_of_edges)
    print("Running time: {:.3f}".format(float(res["time"]))) 


    #original SHD
    #precision, recall, f1-score
    precision = 0  # in the predicted graph, proportion of correct edges
    recall = 0  # among all true edges, proportion of those who are predicted
    n_predicted_edges = 0
    for edge in G_pred.edges:
        if G_pred.has_edge(edge[1], edge[0]):
            if edge[0] < edge[1]:
                n_predicted_edges += 1
        else:
            n_predicted_edges += 1
        if G_true.has_edge(*edge):
            precision += 1
    precision /= n_predicted_edges

    n_true_edges = len(G_true.edges)
    for edge in G_true.edges:
        if G_pred.has_edge(*edge):
            recall += 1
    recall /= n_true_edges

    res["precision"] = precision
    res["recall"] = recall
    res["FDR"] = 1-precision
    res["FNR"] = 1-recall
    res["f1"] = 2 * precision * recall / (precision + recall)
    print("precision, recall, F1 : ", "{:.3f}".format(precision), "{:.3f}".format(recall), "{:.3f}".format(2 * precision * recall / (precision + recall)))



    # ----------------------------------------------------------------------
    # [Option A] Standard MEC F1 (Binary Classification on Ordered Pairs)
    # - Strict: i-j vs i->j counts as partial error (Low Precision)
    # ----------------------------------------------------------------------
    A_true = dag_to_cpdag(G_true_np)
    A_true = (A_true != 0).astype(int)
    np.fill_diagonal(A_true, 0)

    A_pred = (G_pred_np != 0).astype(int)
    np.fill_diagonal(A_pred, 0)

    tp_mec = int(np.logical_and(A_pred == 1, A_true == 1).sum())
    fp_mec = int(np.logical_and(A_pred == 1, A_true == 0).sum())
    fn_mec = int(np.logical_and(A_pred == 0, A_true == 1).sum())

    res["cpdag_precision"] = tp_mec / (tp_mec + fp_mec) if (tp_mec + fp_mec) > 0 else 0.0
    res["cpdag_recall"]    = tp_mec / (tp_mec + fn_mec) if (tp_mec + fn_mec) > 0 else 0.0
    if (res["cpdag_precision"] + res["cpdag_recall"]) > 0:
        res["cpdag_f1"] = 2 * res["cpdag_precision"] * res["cpdag_recall"] / (res["cpdag_precision"] + res["cpdag_recall"])
    else:
        res["cpdag_f1"] = 0.0


    # ----------------------------------------------------------------------
    # [Option B] Benchpress Half-Credit F1
    # - Relaxed: Orientation mismatch gets 0.5 TP / 0.5 FP

    # This is our reported metric in the paper, as F1 score!!
    # ----------------------------------------------------------------------
    def edge_state(A, i, j):
        aij, aji = int(A[i, j] != 0), int(A[j, i] != 0)
        if aij == 0 and aji == 0: return 0  # None
        if aij == 1 and aji == 1: return 1  # Undirected
        if aij == 1 and aji == 0: return 2  # i->j
        return 3                            # j->i

    d = A_true.shape[0]
    TP_b, FP_b, P_b = 0.0, 0.0, 0
    
    # Iterate over unordered pairs {i,j}
    for i in range(d):
        for j in range(i + 1, d):
            st = edge_state(A_true, i, j)
            sp = edge_state(A_pred, i, j)

            if st != 0: P_b += 1   # Count true skeleton edges

            if sp == 0: continue   # Prediction has no edge
            
            if st == 0:            # True has no edge
                FP_b += 1.0        # Full False Positive
            else:                  # Both have edges (Skeleton Match)
                if sp == st:
                    TP_b += 1.0    # Exact Match (Full Credit)
                else: 
                    TP_b += 0.5    # Orientation Mismatch (Half Credit)
                    FP_b += 0.5

    FN_b = float(P_b) - TP_b
    
    # Save as 'benchpress_...' keys
    res["benchpress_precision"] = TP_b / (TP_b + FP_b) if (TP_b + FP_b) > 0 else 0.0
    res["benchpress_recall"]    = TP_b / (TP_b + FN_b) if (TP_b + FN_b) > 0 else 0.0
    
    if (res["benchpress_precision"] + res["benchpress_recall"]) > 0:
        res["benchpress_f1"] = 2 * res["benchpress_precision"] * res["benchpress_recall"] / (res["benchpress_precision"] + res["benchpress_recall"])
    else:
        res["benchpress_f1"] = 0.0

    print("Benchpress F1 (Half-Credit): {:.3f}".format(res["benchpress_f1"]))


    #Ancestor AID


    # ----- Skeleton F1 (Direction Ignored) -----
    skel_pred = G_pred.to_undirected()
    skel_true = G_true.to_undirected()
    
    n_skel_pred = skel_pred.number_of_edges()
    n_skel_true = skel_true.number_of_edges()
    
    tp_skel = 0
    for u, v in skel_pred.edges():
        if skel_true.has_edge(u, v):
            tp_skel += 1
            
    skel_prec = tp_skel / n_skel_pred if n_skel_pred > 0 else 0
    skel_rec = tp_skel / n_skel_true if n_skel_true > 0 else 0
    
    if (skel_prec + skel_rec) > 0:
        skel_f1 = 2 * skel_prec * skel_rec / (skel_prec + skel_rec)
    else:
        skel_f1 = 0.0
        
    res["skeleton_precision"] = skel_prec
    res["skeleton_recall"] = skel_rec
    res["skeleton_f1"] = skel_f1
    print("Skeleton precision, recall, F1 : ", "{:.3f}".format(skel_prec), "{:.3f}".format(skel_rec), "{:.3f}".format(skel_f1))


    out_dir = Path(args.out_root)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _to_scalar(v):
        if isinstance(v, (np.generic,)):
            return v.item()
        return v

    rows = [(k, _to_scalar(v)) for k, v in res.items()]

    df = pd.DataFrame(rows, columns=["metric", "value"])
    df.to_csv(args.stats, index=False, header=False)
    print(f"Saved stats csv -> {args.stats}")

if __name__ == "__main__":
    main()
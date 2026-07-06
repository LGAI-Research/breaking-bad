import argparse
import numpy as np
import pandas as pd
import time
import os
from causallearn.search.ConstraintBased.PC import pc
from causallearn.utils.cit import fisherz, chisq



def convert_causallearn_graph_to_cpdag_adj(cgG):
    """
    causal-learn
      - i -> j : cgG.graph[i,j] = -1 and cgG.graph[j,i] =  1
      - i -- j : cgG.graph[i,j] = -1 and cgG.graph[j,i] = -1
      - i <-> j: cgG.graph[i,j] =  1 and cgG.graph[j,i] =  1  
      - 0 no edge
    """
    Gm = np.asarray(cgG.graph)
    p = Gm.shape[0]
    A = np.zeros((p, p), dtype=int)

    for i in range(p):
        for j in range(i+1, p):
            a = Gm[i, j]
            b = Gm[j, i]
            if a == 0 and b == 0:
                continue

            # undirected i -- j  (-1, -1)
            if a == -1 and b == -1:
                A[i, j] = 1
                A[j, i] = 1

            # i -> j  (-1, +1) on (i,j) / (j,i)
            elif a == -1 and b == 1:
                A[i, j] = 1
                A[j, i] = 0

            # j -> i  (+1, -1)
            elif a == 1 and b == -1:
                A[j, i] = 1
                A[i, j] = 0

            # bidirected i <-> j  (+1, +1)
            elif a == 1 and b == 1:

                A[i, j] = 0
                A[j, i] = 0
                print(f"[warn] bidirected edge {i}<->{j} detected (PC usually shouldn't output this).")

            else:
                A[i, j] = 0
                A[j, i] = 0
                # print(f"[warn] unexpected encoding ({a},{b}) for edge {i}-{j}; treating as undirected.")

    np.fill_diagonal(A, 0)
    return A


def main():
    parser = argparse.ArgumentParser(description="Run PC Algorithm (causal-learn)")
    parser.add_argument("--input", type=str, required=True, help="Path to input .npy data file")
    parser.add_argument("--output", type=str, required=True, help="Path to save output result (.csv)")
    parser.add_argument("--stats", type=str, required=False, help="Path to save stats (.csv)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level (default 0.05)")
    args = parser.parse_args()


    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    data = np.load(args.input)
    
    start_time = time.time()
    

    # cg = pc(data, args.alpha, fisherz, True, 0, -1) 
    cg = pc(data, alpha=args.alpha, indep_test="fisherz", stable=True, uc_rule=0, uc_priority=2, show_progress=True)
    
    
    Gm = np.asarray(cg.G.graph)
    vals, cnts = np.unique(Gm, return_counts=True)
    print("cg.G.graph unique values:", dict(zip(vals, cnts)))


    end_time = time.time()
    elapsed = end_time - start_time


    adj_mat = convert_causallearn_graph_to_cpdag_adj(cg.G)
    
    df_res = pd.DataFrame(adj_mat)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df_res.to_csv(args.output, index=False)
    
    if args.stats:
        with open(args.stats, 'w') as f:
  
            f.write(f"time,{elapsed:.6f}\n")
            f.write(f"algorithm,PC\n")
            f.write(f"nodes,{data.shape[1]}\n")
            f.write(f"alpha,{args.alpha}\n")
    
    print(f"[PC] Completed. Time: {elapsed:.4f}s, Nodes: {data.shape[1]}, Alpha: {args.alpha}")

if __name__ == "__main__":
    main()
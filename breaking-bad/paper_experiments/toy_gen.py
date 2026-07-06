import numpy as np
import random
from sempler.generators import dag_avg_deg
from utils_lges import get_random_data, dag_to_cpdag

def main():
    # seed
    np.random.seed(451)
    random.seed(451)
    G = dag_avg_deg(5, 2.5, 1, 1, random_state=451)
    
    # CPDAG
    cpdag = dag_to_cpdag(G)
    
    np.random.seed(451)
    random.seed(451)
    X = get_random_data(G, 2000, model='linear-gaussian', interventions=[{}], random_state=451)[0]
    
    np.save('toy_data.npy', np.ascontiguousarray(X.astype(np.float64)))
    np.savetxt('toy_gt_dag.csv', G, delimiter=',', fmt='%d', comments='')
    np.savetxt('toy_gt.csv', cpdag, delimiter=',', fmt='%d', comments='') 
    
    print("✅ toy_data.npy generated.")
    print("\n[Ground Truth DAG]")
    print(G)
    print("\n[Ground Truth CPDAG]")
    print(cpdag)

if __name__ == "__main__":
    main()
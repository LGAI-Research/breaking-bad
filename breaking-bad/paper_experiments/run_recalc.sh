# ==============================
# Parameter Lists for repetitive experiments
# ==============================
#VAR_LIST=(10 20 50 100 200 300 400 500 700 1000)   
VAR_LIST=(20)              # ex) number of variables
DEG_LIST=(2.0)              # ex) average degree for ER graphs
#DEG_LIST=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0)      
SAMPLES_LIST=(1000)      # ex) number of samples
#SAMPLES_LIST=(100 200 500 1000 10000 100000 1000000)
 
VARIANTS_LIST=(0 1 2 3) # 0: GES, 1: GES-D, 2: GES-DP, 3: GES-DP+
#ALPHA_LIST=(1 2 3 4 7 10 12 15 18 20 30)    # ex) BIC hyperparameter alpha
#ALPHA_LIST=(0.05)
ALPHA_LIST=(2)


#BASELINE_LIST=("ges" "ops" "xges")           # "ges" "ops" "xges" "boss" "pc"
#BASELINE_LIST=("ops" "xges" "boss" "lges-safe" "lges-cons") \\
BASELINE_LIST=("ges")   


M_LIST=(2) #scale-free density parameter
# ==============================
# fixed parameters
# ==============================
DELETE_OP=1
RUNS=30

GRAPH_TYPE="ER" # Erdos-Renyi
#GRAPH_TYPE="SF" # Scale-Free

echo "=================================================="
echo "🔄 Starting Stats Re-calculation (Checking Skeleton F1)"
echo "=================================================="

for VAR in "${VAR_LIST[@]}"; do
  for DEG in "${DEG_LIST[@]}"; do
    for samples in "${SAMPLES_LIST[@]}"; do
      for alpha in "${ALPHA_LIST[@]}"; do
        for variant in "${VARIANTS_LIST[@]}"; do
          for baseline in "${BASELINE_LIST[@]}"; do

            # python recalc_stats.py
            python recalc_stats.py \
              --var "$VAR" \
              --deg "$DEG" \
              --samples "$samples" \
              --variant "$variant" \
              --delete_op "$DELETE_OP" \
              --alpha "$alpha" \
              --runs "$RUNS" \
              --baseline "$baseline" \
              --graph_type "$GRAPH_TYPE"

          done
        done
      done
    done
  done
done

echo ""
echo "=================================================="
echo "🎉 Stats updated successfully!"
echo "=================================================="
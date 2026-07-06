# ==============================
# Parameter Lists for repetitive experiments
# ==============================
#VAR_LIST=(10 20 50 100 200 300 400 500)   
VAR_LIST=(100)              # ex) number of variables
DEG_LIST=(5.0)              # ex) average degree for ER graphs
#DEG_LIST=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0)
SAMPLES_LIST=(10000)      # ex) number of samples
#SAMPLES_LIST=(1000 10000 100000)
#SAMPLES_LIST=(100 200 500 1000 10000 100000 1000000)

VARIANTS_LIST=(0 1 2 3)  
#VARIANTS_LIST=(0) # 0: GES, 1: GES-D, 2: GES-DP, 3: GES-DP+, 4: GES-DP0
ALPHA_LIST=(1 2 3 4 7 10 12 15 18 20 30)    # ex) BIC hyperparameter alpha
#ALPHA_LIST=(2)
#ALPHA_LIST=(1 4)

BASELINE_LIST=("ges" "xges" "ops" "boss" "lges-safe" "lges-cons")
#BASELINE_LIST=("boss")   


M_LIST=(2) #scale-free density parameter
# ==============================
# fixed parameters
# ==============================
DELETE_OP=1
RUNS=30

GRAPH_TYPE="ER" # Erdos-Renyi
#GRAPH_TYPE="SF" # Scale-Free
#GRAPH_TYPE="NL_ER" # Non-linear ER

EVAL_ONLY=1 # Set to 1 to only evaluate existing results without running new experiments


echo "=================================================="
echo "🚀 Starting Large Scale Batch Experiment"
echo "=================================================="

# ==============================
# Sextuple Loop
# Var → Deg → Samples → Alpha → Variant → Baseline
# ==============================
for VAR in "${VAR_LIST[@]}"; do
  for DEG in "${DEG_LIST[@]}"; do
    for samples in "${SAMPLES_LIST[@]}"; do
      for alpha in "${ALPHA_LIST[@]}"; do
        for variant in "${VARIANTS_LIST[@]}"; do
          for baseline in "${BASELINE_LIST[@]}"; do
            for m in "${M_LIST[@]}"; do

              echo ""
              echo "--------------------------------------------------"
              echo "▶ Running Configuration:"
              echo "   Variables : $VAR"
              echo "   Degree    : $DEG"
              echo "   Samples   : $samples"
              echo "   Alpha     : $alpha"
              echo "   Variant   : $variant"
              echo "   Baseline  : $baseline"
              echo "   M         : $m"
              echo "--------------------------------------------------"

              EVAL_FLAG=""
              if [ "$EVAL_ONLY" -eq 1 ]; then
                EVAL_FLAG="--eval_only"
              fi

              python batch_run.py \
                --var "$VAR" \
                --deg "$DEG" \
                --samples "$samples" \
                --variant "$variant" \
                --delete_op "$DELETE_OP" \
                --alpha "$alpha" \
                --runs "$RUNS" \
                --baseline "$baseline" \
                --graph_type "$GRAPH_TYPE"\
                --m "$m"\
                $EVAL_FLAG

              if [ $? -eq 0 ]; then
                echo "✅ Completed (V=$VAR, D=$DEG, S=$samples, A=$alpha, Vt=$variant, B=$baseline)"
              else
                echo "❌ Failed    (V=$VAR, D=$DEG, S=$samples, A=$alpha, Vt=$variant, B=$baseline)"
              fi
            done
          done
        done
      done
    done
  done
done

echo ""
echo "=================================================="
echo "🎉 All experiments finished!"
echo "=================================================="

DATASET_LIST=("ARTH150")     # "MAGIC-IRRI" "Sachs"
SAMPLES_LIST=(10000)            
ALPHA_LIST=(2.0)                # Score Alpha
#ALPHA_LIST=(0.05)  
VARIANTS_LIST=(0 1 2 3)         # 0:Base, 1:Single, 2:DP, 3:DP+
BASELINE_LIST=("ges")   # "ops" "xges" "lges-cons" "lges-safe" "boss"


RUNS=30  
VERBOSE=1

echo "=================================================="
echo "🚀 Starting Real-world Batch Experiment"
echo "=================================================="

# ==============================
# Quadruple Loop
# Dataset → Samples → Alpha → Variant → Baseline
# ==============================
for DATASET in "${DATASET_LIST[@]}"; do
    for samples in "${SAMPLES_LIST[@]}"; do
      for alpha in "${ALPHA_LIST[@]}"; do
        for variant in "${VARIANTS_LIST[@]}"; do
          for baseline in "${BASELINE_LIST[@]}"; do

            echo ""
            echo "--------------------------------------------------"
            echo "▶ Running Configuration:"
            echo "   Dataset   : $DATASET"
            echo "   Samples   : $samples"
            echo "   Alpha     : $alpha"
            echo "   Variant   : $variant"
            echo "   Baseline  : $baseline"
            echo "--------------------------------------------------"

            python batch_run_realworld.py \
              --dataset "$DATASET" \
              --samples "$samples" \
              --variant "$variant" \
              --alpha "$alpha" \
              --runs "$RUNS" \
              --baseline "$baseline" \
              --verbose "$VERBOSE"

            if [ $? -eq 0 ]; then
              echo "✅ Completed (Data=$DATASET, S=$samples, A=$alpha, Vt=$variant, B=$baseline)"
            else
              echo "❌ Failed    (Data=$DATASET, S=$samples, A=$alpha, Vt=$variant, B=$baseline)"
            fi

          done
        done
      done
    done
done

echo ""
echo "=================================================="
echo "🎉 All Real-world experiments finished!"
echo "=================================================="
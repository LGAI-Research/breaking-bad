#!/bin/bash

# ==============================
# Parameter Lists
# ==============================
DATASET_LIST=("MAGIC-NIAB")
SAMPLES_LIST=(10000)
#SAMPLES_LIST=(853)
ALPHAS=(2.0)
VARIANTS=(0 1 2 3)
#BASELINES=("ges" "xges" "ops" "boss" "lges-safe" "lges-cons") 
BASELINES=("boss") 
RUNS=30
#RUNS=1
VERBOSE=1

echo "=================================================="
echo "🚀 Starting Realworld Batch Experiment"
echo "=================================================="

fail_count=0

# ==============================
# Nested Loop
# Dataset → Samples → Alpha → Variant → Baseline
# ==============================
for dataset in "${DATASET_LIST[@]}"; do
    for samples in "${SAMPLES_LIST[@]}"; do
        for baseline in "${BASELINES[@]}"; do
            for alpha in "${ALPHAS[@]}"; do
                for variant in "${VARIANTS[@]}"; do
                    
                    echo ""
                    echo "--------------------------------------------------------"
                    echo "▶️ Running Experiment:"
                    echo "   Dataset   : $dataset"
                    echo "   Samples   : $samples"
                    echo "   Baseline  : $baseline"
                    echo "   Variant   : $variant"
                    echo "   Alpha     : $alpha"
                    echo "--------------------------------------------------------"

                    python batch_run_realworld.py \
                        --dataset "$dataset" \
                        --samples "$samples" \
                        --alpha "$alpha" \
                        --variant "$variant" \
                        --baseline "$baseline" \
                        --runs "$RUNS" \
                        --verbose "$VERBOSE"

                    if [ $? -eq 0 ]; then
                        echo "✅ Success: $dataset | $baseline | V=$variant | A=$alpha"
                    else
                        echo "❌ Failed : $dataset | $baseline | V=$variant | A=$alpha"
                        ((fail_count++))
                    fi
                    
                done
            done
        done
    done
done

echo "========================================================"
if [ $fail_count -gt 0 ]; then
    echo "⚠️ Experiments finished with $fail_count failures."
else
    echo "🎉 All experiments completed successfully!"
fi
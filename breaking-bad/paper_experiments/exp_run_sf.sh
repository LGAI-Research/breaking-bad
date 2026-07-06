
mkdir -p ./outputs

var=50
samples=1000
seed=3      
verbose_level=2
variant_type=1
delete_op=1
alpha=2.0
baseline_type="ges"
viz_log=false

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --var) var="$2"; shift ;;
        --samples) samples="$2"; shift ;;
        --seed) seed="$2"; shift ;;        # graph_num -> seed
        --alpha) alpha="$2"; shift ;;
        --baseline) baseline_type="$2"; shift ;; 
        --variant) variant_type="$2"; shift ;;
        --verbose) verbose_level="$2"; shift ;;
        *) echo "Unknown parameter: $1"; shift ;; 
    esac
    shift
done


# GT: ScaleFree_DAG_var{num_var}_seed{seed}.csv
ground_truth="ScaleFree_DAG_var${var}_seed${seed}.csv"

# Data: ScaleFree_Data_var{num_var}_n{num_samples}_seed{seed}.npy
data="ScaleFree_Data_var${var}_n${samples}_seed${seed}.npy"

# Tag : run-sf-var50-n1000-seed1-...
tag="run-sf-var${var}-n${samples}-seed${seed}-alpha${alpha}-${baseline_type}-var${variant_type}"


RESULT_DIR="./outputs/scalefree/${baseline_type}/var${var}"
mkdir -p "$RESULT_DIR"

VIZ_OPT=""
# if [ "$viz_log" = true ]; then ...

BASELINE_OPT=""
if [ "$baseline_type" != "" ]; then
    BASELINE_OPT="--baseline $baseline_type"
fi

echo "▶ Running Scale-Free Exp: var=$var, seed=$seed, algo=$baseline_type"

# ------------------------------------------------------------------
# (C++)
# ------------------------------------------------------------------
../build/src/Search \
    --input "../data/samples/$data" \
    -a ${alpha} \
    --output "${RESULT_DIR}/${tag}_result.csv" \
    --stats "${RESULT_DIR}/${tag}_stats.csv" \
    -v ${verbose_level} \
    --variant ${variant_type} \
    --delete_op ${delete_op} \
    $BASELINE_OPT

if [ $? -ne 0 ]; then
    echo "❌ Error: execution failed for $tag"
    exit 1
fi

# ------------------------------------------------------------------
#  (Python)
# ------------------------------------------------------------------
python simple_exp.py \
    --result "${RESULT_DIR}/${tag}_result.csv" \
    --stats "${RESULT_DIR}/${tag}_stats.csv" \
    --ground_truth "../data/ground_truth/$ground_truth" \
    --tag "$tag" \
    --out-root "$RESULT_DIR"

echo "✅ Done: $tag"
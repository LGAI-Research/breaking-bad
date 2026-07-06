var=20
deg=2.0
samples=1000
verbose_level=2
variant_type=0
delete_op=1
graph_num=30
alpha=2.0
baseline_type="ges"
graph_type="ER"  
m_value=2

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --var) var="$2"; shift ;;
        --deg) deg="$2"; shift ;;
        --samples) samples="$2"; shift ;;
        --verbose) verbose_level="$2"; shift ;;
        --variant) variant_type="$2"; shift ;;
        --delete_op) delete_op="$2"; shift;;
        --baseline) baseline_type="$2"; shift ;; 
        --alpha) alpha="$2"; shift ;;            
        --graph_num) graph_num="$2"; shift ;;
        --graph_type) graph_type="$2"; shift ;;    
        --m) m_value="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$graph_type" ]; then
    graph_type="ER"
fi

if [ "$graph_type" == "SF" ]; then
    ground_truth="ScaleFree_m${m_value}_DAG_var${var}_seed${graph_num}.csv"
    data="ScaleFree_m${m_value}_Data_var${var}_n${samples}_seed${graph_num}.npy"
#non-linear ER graph
elif [ "$graph_type" == "NL_ER" ]; then
    ground_truth="True_DAG_nonlinear_var${var}_avg_deg${deg}_graph_num${graph_num}.csv"
    data="Data_nonlinear_var${var}_avg_deg${deg}_n_samples${samples}_graph_num${graph_num}.npy"
# Gumbel ER graph
elif [ "$graph_type" == "gumbel" ]; then
    ground_truth="True_DAG_gumbel_var${var}_avg_deg${deg}_graph_num${graph_num}.csv"
    data="Data_gumbel_var${var}_avg_deg${deg}_n_samples${samples}_graph_num${graph_num}.npy"
else
    ground_truth="True_DAG_var${var}_avg_deg${deg}_graph_num${graph_num}.csv"
    data="Data_var${var}_avg_deg${deg}_n_samples${samples}_graph_num${graph_num}.npy"
fi

if [ "$graph_type" == "SF" ]; then
    tag=run-${graph_type}-m${m_value}-${var}-${deg}-${samples}-alpha${alpha}-${baseline_type}-variant${variant_type}-graph${graph_num}-search
else
    tag=run-${graph_type}-${var}-${deg}-${samples}-alpha${alpha}-${baseline_type}-variant${variant_type}-graph${graph_num}-search
fi


BASELINE_OPT=""
if [ -n "$baseline_type" ]; then
    BASELINE_OPT="--baseline $baseline_type"
fi

echo "Running with: var=$var, alpha=$alpha, variant=$variant_type, baseline=$baseline_type, graph=$graph_num"


if [ "$graph_type" == "SF" ]; then
    RESULT_DIR="./outputs/${baseline_type}/${graph_type}/m${m_value}/var${var}/deg${deg}"
else
    RESULT_DIR="./outputs/${baseline_type}/${graph_type}/var${var}/deg${deg}"
fi

# Ensure the result directory exists before the search writes into it
mkdir -p "${RESULT_DIR}"



# ------------------------------------------------------------------
# C++
# ------------------------------------------------------------------
../build/src/Search --input ../data/samples/$data \
       -a ${alpha} \
       --output ${RESULT_DIR}/${tag}_result.csv \
       --stats ${RESULT_DIR}/${tag}_stats.csv \
       -v ${verbose_level} \
       --variant ${variant_type} \
       --delete_op ${delete_op} \
       $BASELINE_OPT

if [ $? -ne 0 ]; then
    echo "❌ Error: XGES execution failed."
    exit 1
fi 

# ------------------------------------------------------------------
# Python evaluation
# ------------------------------------------------------------------
python simple_exp.py \
       --result ${RESULT_DIR}/${tag}_result.csv \
       --stats ${RESULT_DIR}/${tag}_stats.csv \
       --ground_truth ../data/ground_truth/${ground_truth} \
       --tag $tag \
       --out-root ./outputs
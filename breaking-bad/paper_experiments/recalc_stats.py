import subprocess
import pandas as pd
import argparse
import os


def parse_ecshd(val):
    if isinstance(val, str):
        clean_val = val.strip("()")
        try:
            return float(clean_val.split(",")[0])
        except:
            return None
    elif isinstance(val, (int, float)):
        return float(val)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", type=int, default=100)
    parser.add_argument("--deg", type=float, default=2.0)
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--variant", type=int, default=2)
    parser.add_argument("--delete_op", type=int, default=1)
    parser.add_argument("--alpha", type=float, default=2.0)
    parser.add_argument("--runs", type=int, default=30)  
    parser.add_argument("--baseline", type=str, default="ges")
    parser.add_argument("--graph_type", type=str, default="ER")
    parser.add_argument("--m", type=int, default=4) 

    args = parser.parse_args()

    results = []

    print(f"🔄 Re-calculating Stats: type={args.graph_type}, var={args.var}, runs={args.runs}")

    for i in range(1, args.runs + 1):

        if args.graph_type == "SF":
            output_dir = f"./outputs/{args.baseline}/{args.graph_type}/m{args.m}/var{args.var}/deg{args.deg}"
            tag = f"run-{args.graph_type}-m{args.m}-{args.var}-{args.deg}-{args.samples}-alpha{args.alpha}-{args.baseline}-variant{args.variant}-graph{i}-search"
            gt_file = f"../data/ground_truth/ScaleFree_DAG_var{args.var}_seed{i}.csv"
        else:
            output_dir = f"./outputs/{args.baseline}/{args.graph_type}/var{args.var}/deg{args.deg}"
            tag = f"run-{args.graph_type}-{args.var}-{args.deg}-{args.samples}-alpha{args.alpha}-{args.baseline}-variant{args.variant}-graph{i}-search"
            gt_file = f"../data/ground_truth/True_DAG_var{args.var}_avg_deg{args.deg}_graph_num{i}.csv"

        output_file = f"{output_dir}/{tag}_result.csv"
        stats_file = f"{output_dir}/{tag}_stats.csv"

        if not os.path.exists(output_file):
            print(f"⚠️ Missing result file: {output_file}")
            continue

        # simple_exp
        eval_cmd = [
            "python", "simple_exp.py",
            "--result", output_file,
            "--stats", stats_file,
            "--ground_truth", gt_file,
            "--tag", tag,
            "--out-root", "./outputs"
        ]

        try:
            subprocess.run(eval_cmd, check=True, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"❌ Error computing stats for run {i}")
            continue

        # stats
        if os.path.exists(stats_file):
            try:
                df = pd.read_csv(stats_file, header=None, names=["metric", "value"])
                keys = df["metric"].astype(str).str.strip()
                res_dict = dict(zip(keys, df["value"]))
                res_dict["graph_num"] = i
                results.append(res_dict)
            except Exception as e:
                print(f"⚠️ Error reading {stats_file}: {e}")

    # ============================
    # Summary 
    # ============================
    if not results:
        print("❌ No results collected.")
        return

    df_res = pd.DataFrame(results)


    df_res["graph_num"] = pd.to_numeric(df_res["graph_num"], errors="coerce")
    df_res = df_res.sort_values("graph_num").head(args.runs)

    if "ecshd" in df_res.columns:
        df_res["ecshd"] = df_res["ecshd"].apply(parse_ecshd)

    for c in df_res.columns:
        df_res[c] = pd.to_numeric(df_res[c], errors="coerce")
    
    df_res["graph_num"] = pd.to_numeric(df_res["graph_num"], errors="coerce")
    df_res = df_res.dropna(subset=["graph_num"])

    df_res = df_res.sort_values("graph_num")

    df_res = df_res[df_res["graph_num"] <= 30]

    if len(df_res) > 30:
        df_res = df_res.head(30)

    numeric_cols = df_res.select_dtypes(include=["number"]).columns.tolist()
    summary = df_res[numeric_cols].describe().T[["mean", "std"]]

    if args.graph_type == "SF":
        summary_dir = f"./outputs/summary/{args.baseline}/{args.graph_type}/m{args.m}/var{args.var}/deg{args.deg}"
        filename = f"summary_{args.graph_type}_m{args.m}_var{args.var}_deg{args.deg}_samples{args.samples}_alpha{args.alpha}_{args.baseline}_variant{args.variant}.csv"
    else:
        summary_dir = f"./outputs/summary/{args.baseline}/{args.graph_type}/var{args.var}/deg{args.deg}"
        filename = f"summary_{args.graph_type}_var{args.var}_deg{args.deg}_samples{args.samples}_alpha{args.alpha}_{args.baseline}_variant{args.variant}.csv"

    os.makedirs(summary_dir, exist_ok=True)
    save_path = os.path.join(summary_dir, filename)

    df_res.to_csv(save_path, index=False)

    with open(save_path, "a") as f:
        f.write("\nSummary (Mean/Std)\n")
        summary.to_csv(f)

    print("\n📈 Batch Summary (Mean ± Std)")
    print(summary)
    print(f"\n✅ Saved to: {save_path}")


if __name__ == "__main__":
    main()
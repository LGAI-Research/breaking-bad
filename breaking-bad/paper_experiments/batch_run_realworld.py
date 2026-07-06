import subprocess
import pandas as pd
import argparse
import os
import sys

def parse_ecshd(val):
    if isinstance(val, str):
        clean_val = val.strip("()")
        try:
            parts = clean_val.split(",")
            return float(parts[0])
        except (ValueError, IndexError):
            return None
    elif isinstance(val, (int, float)):
        return float(val)
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, default="MAGIC-IRRI")
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--variant", type=int, default=2)
    parser.add_argument("--alpha", type=float, default=2.0)
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--baseline", type=str, default="ges")
    parser.add_argument("--verbose", type=int, default=1)
    
    args = parser.parse_args()

    INPUT_ROOT = "../data/Realworld"
    BUILD_PATH = "../build/src/Search"
    
    results = []

    print(f"🚀 Starting Real-world Batch: dataset={args.dataset}, samples={args.samples}, alpha={args.alpha}, variant={args.variant}, baseline='{args.baseline}'")

    for i in range(1, args.runs + 1):
        
        # Output Directory: ./outputs/Realworld/{dataset}/{baseline}
        output_dir = f"./outputs/Realworld/{args.dataset}/{args.baseline}"
        os.makedirs(output_dir, exist_ok=True)
        
        tag = f"run-{args.dataset}_{args.samples}-seed{i}-alpha{args.alpha}-{args.baseline}-variant{args.variant}-search"
        
        # File Paths
        input_file = f"{INPUT_ROOT}/{args.dataset}/Data_{args.dataset}_{args.samples}_{i}.npy"
        gt_file = f"{INPUT_ROOT}/{args.dataset}/True_DAG_{args.dataset}.csv"
        
        output_file = f"{output_dir}/{tag}_result.csv"
        stats_file = f"{output_dir}/{tag}_stats.csv"

        if not os.path.exists(input_file):

            input_file_lower = f"{INPUT_ROOT}/{args.dataset}/data_{args.dataset}_{args.samples}_{i}.npy"
            if os.path.exists(input_file_lower):
                input_file = input_file_lower
            else:
                print(f"⚠️ Input file not found: {input_file}. Skipping run {i}.")
                continue

        # ==========================================
        # Baseline (PC vs XGES/Others)
        # ==========================================
        if args.baseline == "pc":
            pc_alpha = 0.05 if args.alpha >= 1.0 else args.alpha
            cmd = [
                "python", "run_pc.py",
                "--input", input_file,
                "--output", output_file,
                "--stats", stats_file,
                "--alpha", str(pc_alpha)
            ]
        else:
            # XGES / GES / OPS / BOSS 
            baseline_arg = []
            if args.baseline:
                 baseline_arg = ["--baseline", args.baseline]

            cmd = [
                BUILD_PATH,
                "--input", input_file,
                "--output", output_file,
                "--stats", stats_file,
                "-a", str(args.alpha),
                "--variant", str(args.variant),
                "-v", str(args.verbose)
            ] + baseline_arg

            # run-{dataset}_{samples}-seed{run_id}-alpha{alpha}-boss-variant0-search_result.csv
            if args.baseline == "boss":
                boss_msg = ""
                # Variant 0 Result File Path Construction
                v0_filename = f"run-{args.dataset}_{args.samples}-seed{i}-alpha{args.alpha}-boss-variant0-search_result.csv"
                v0_path = os.path.join(output_dir, v0_filename)
 
                if os.path.exists(v0_path):
                    cmd += ["--initial_pdag", v0_path]
                    boss_msg = f"(Found V0: {v0_filename})"
                else:
             
                    print(f"⚠️ BOSS Variant {args.variant} requires Variant 0 result, but not found: {v0_path}")
    

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            print(f"❌ Error occurred at run {i} (Cmd: {cmd[0]})")
            continue

        eval_cmd = [
            "python", "simple_exp.py",
            "--result", output_file,
            "--stats", stats_file,
            "--ground_truth", gt_file,
            "--tag", tag,
            "--out-root", "./outputs"
        ]
        
        try:
             subprocess.run(eval_cmd, check=True)
        except subprocess.CalledProcessError:
             pass

        if os.path.exists(stats_file):
            try:
                df = pd.read_csv(stats_file, header=None, names=["metric", "value"])
                keys = df["metric"].astype(str).str.strip()
                res_dict = dict(zip(keys, df["value"]))
                res_dict["graph_num"] = i # Seed
                results.append(res_dict)
                print(f"   [Run {i}] SHD: {res_dict.get('shd')}, F1: {res_dict.get('f1')}")
            except Exception as e:
                print(f"⚠️ Error reading Stats CSV {stats_file}: {e}")
        else:
            print(f"⚠️ Stats file not found: {stats_file}")

    if results:
        df_res = pd.DataFrame(results)
        
        if "ecshd" in df_res.columns:
            df_res["ecshd"] = df_res["ecshd"].apply(parse_ecshd)

        for c in df_res.columns:
            try:
                df_res[c] = pd.to_numeric(df_res[c], errors='coerce') 
            except ValueError:
                pass
        
        summary_dir = f"./outputs/summary/Realworld/{args.dataset}"
        os.makedirs(summary_dir, exist_ok=True)
        
        save_filename = f"summary_{args.dataset}_samples{args.samples}_alpha{args.alpha}_{args.baseline}_variant{args.variant}.csv"
        save_path = os.path.join(summary_dir, save_filename)
        
        print_cols = ["time", "shd", "ecshd", "f1", "precision", "recall", "benchpress_f1"]
        existing_print_cols = [c for c in print_cols if c in df_res.columns]
        
        summary = df_res.describe().T[["mean", "std"]]
        
        print("\n" + "-"*30)
        print(f"📈 Batch Summary: {args.dataset} (Mean ± Std)")
        if existing_print_cols:
             print(summary.loc[existing_print_cols])
        print("-"*30)

        expected_cols = ["graph_num"] + existing_print_cols
        remaining_cols = [c for c in df_res.columns if c not in expected_cols]
        final_cols = expected_cols + remaining_cols
        
        df_res = df_res.reindex(columns=final_cols)
        df_res.to_csv(save_path, index=False)
        
        with open(save_path, 'a') as f:
            f.write("\nSummary (Mean/Std)\n")
            summary.to_csv(f)

        print(f"✅ Summary saved to: {save_path}")

    else:
        print("\n❌ No results collected properly.")

if __name__ == "__main__":
    main()
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
            return float(parts[0]) # (SHD)
        except (ValueError, IndexError):
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
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--baseline", type=str, default="ges")
    #parser.add_argument("--graph_type", type=str, default="ER", help="Graph type: 'ER' or 'SF'")
    parser.add_argument("--graph_type", type=str, default="ER", help="Graph type: 'ER', 'SF', or 'NL_ER'")
    parser.add_argument("--m", type=int, default=4) # 
    parser.add_argument("--eval_only", action="store_true", help="Skip graph discovery and only evaluate existing results")
    
    args = parser.parse_args()

    results = []

    print(f"🚀 Starting Batch: type={args.graph_type}, var={args.var}, samples={args.samples}, alpha={args.alpha}, variant={args.variant}, baseline='{args.baseline}'")

    for i in range(1, args.runs + 1):
        

        if args.graph_type == "SF":
            output_dir = f"./outputs/{args.baseline}/{args.graph_type}/m{args.m}/var{args.var}/deg{args.deg}"
        # Non-linear ER graph
        else:
            output_dir = f"./outputs/{args.baseline}/{args.graph_type}/var{args.var}/deg{args.deg}"
        os.makedirs(output_dir, exist_ok=True)

        if args.graph_type == "SF":
            tag = f"run-{args.graph_type}-m{args.m}-{args.var}-{args.deg}-{args.samples}-alpha{args.alpha}-{args.baseline}-variant{args.variant}-graph{i}-search"
            gt_file = f"../data/ground_truth/ScaleFree_m{args.m}_DAG_var{args.var}_seed{i}.csv"
        else:
            tag = f"run-{args.graph_type}-{args.var}-{args.deg}-{args.samples}-alpha{args.alpha}-{args.baseline}-variant{args.variant}-graph{i}-search"
            gt_file = f"../data/ground_truth/True_DAG_var{args.var}_avg_deg{args.deg}_graph_num{i}.csv"
        
        output_file = f"{output_dir}/{tag}_result.csv"
        stats_file = f"{output_dir}/{tag}_stats.csv"

        if not args.eval_only:
            if args.baseline == "pc":
                if args.graph_type == "SF":
                    input_file = f"../data/samples/ScaleFree_m{args.m}_Data_var{args.var}_n{args.samples}_seed{i}.npy"
                else:
                    input_file = f"../data/samples/Data_var{args.var}_avg_deg{args.deg}_n_samples{args.samples}_graph_num{i}.npy"
                
                pc_alpha = 0.05 if args.alpha >= 1.0 else args.alpha
                
                cmd = [
                    "python", "run_pc.py",
                    "--input", input_file,
                    "--output", output_file,
                    "--stats", stats_file,
                    "--alpha", str(pc_alpha),
                ]
                
                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError:
                    print(f"❌ Error occurred at run {i} (Cmd: {cmd[0]})")
                    continue
                    
            else:
                cmd = [
                    "bash", "exp_run.sh",
                    "--var", str(args.var),
                    "--deg", str(args.deg),
                    "--samples", str(args.samples),
                    "--variant", str(args.variant),
                    "--delete_op", str(args.delete_op),
                    "--alpha", str(args.alpha),
                    "--graph_num", str(i),
                    "--verbose", "1",
                    "--baseline", args.baseline,
                    "--graph_type", args.graph_type,   
                    "--m", str(args.m)  # scale-free
                ]

                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError:
                    print(f"❌ Error occurred at run {i} (Cmd: {cmd[0]})")
                    continue

        # 강제로 모든 알고리즘에 대해 (eval_only 모드이거나 pc 모드일 경우 등) 평가 실행
        if args.eval_only or args.baseline == "pc":
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
            except Exception as e:
                 print(f"⚠️ Error in evaluation: {e}")

        stats_path = f"{output_dir}/{tag}_stats.csv"

        # if args.baseline == "pc":
        #     if args.graph_type == "SF":
        #          gt_file = f"../data/ground_truth/ScaleFree_DAG_var{args.var}_seed{i}.csv"
        #     else:
        #          gt_file = f"../data/ground_truth/True_DAG_var{args.var}_avg_deg{args.deg}_graph_num{i}.csv"
            
        #     eval_cmd = [
        #         "python", "simple_exp.py",
        #         "--result", output_file,
        #         "--stats", stats_file,
        #         "--ground_truth", gt_file,
        #         "--tag", tag,
        #         "--out-root", "./outputs"
        #     ]
        #     try:
        #          subprocess.run(eval_cmd, check=True)
        #     except:
        #          pass


        stats_path = f"{output_dir}/{tag}_stats.csv"

        if os.path.exists(stats_path):
            try:
                df = pd.read_csv(stats_path, header=None, names=["metric", "value"])
                keys = df["metric"].astype(str).str.strip()
                res_dict = dict(zip(keys, df["value"]))
                res_dict["graph_num"] = i
                results.append(res_dict)
                print(f"   [Run {i}] SHD: {res_dict.get('shd')}, F1: {res_dict.get('f1')}")
            except Exception as e:
                print(f"⚠️ Error reading CSV {stats_path}: {e}")
        else:
            print(f"⚠️ Stats file not found: {stats_path}")


    if results:
        df_res = pd.DataFrame(results)
        
        if "ecshd" in df_res.columns:
            df_res["ecshd"] = df_res["ecshd"].apply(parse_ecshd)

        expected_cols = [
            "graph_num", "time", 
            "shd", "ecshd",
            "order_divergence", 
            "precision", "recall", "f1", 
            "cpdag_precision", "cpdag_recall", "cpdag_f1", 
            "skeleton_precision", "skeleton_recall", "skeleton_f1",
            "benchpress_precision", "benchpress_recall", "benchpress_f1"
            "ancestor_aid_norm", "ancestor_aid_mistakes"
        ]
        

        existing_cols = [c for c in expected_cols if c in df_res.columns]
        other_cols = [c for c in df_res.columns if c not in expected_cols]
        final_cols = existing_cols + other_cols 


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
        else:
            summary_dir = f"./outputs/summary/{args.baseline}/{args.graph_type}/var{args.var}/deg{args.deg}"

        os.makedirs(summary_dir, exist_ok=True)


        if args.graph_type == "SF":
            save_filename = (
                f"summary_{args.graph_type}_m{args.m}_var{args.var}_deg{args.deg}_"
                f"samples{args.samples}_alpha{args.alpha}_{args.baseline}_variant{args.variant}.csv"
            )
        else:
            save_filename = (
                f"summary_{args.graph_type}_var{args.var}_deg{args.deg}_"
                f"samples{args.samples}_alpha{args.alpha}_{args.baseline}_variant{args.variant}.csv"
            )

        save_path = os.path.join(summary_dir, save_filename)
        

        print_cols = [c for c in existing_cols if c in summary.index and c != "graph_num"]

        print("\n" + "-"*30)
        print("📈 Batch Summary (Mean ± Std)")
        print(summary.loc[print_cols]) 
        print("-"*30)

        df_res = df_res.reindex(columns=final_cols)
        df_res.to_csv(save_path, index=False)
        
        with open(save_path, 'a') as f:
            f.write("\nSummary (Mean/Std)\n")
            summary.loc[print_cols].to_csv(f)

        print(f"✅ Summary saved to: {save_path}")
    else:
        print("\n❌ No results collected properly.")

if __name__ == "__main__":
    main()
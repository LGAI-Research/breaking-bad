import os
import re
import pandas as pd
import numpy as np
from io import StringIO


# =========================================================
# 1. LOAD SUMMARY STATS
# =========================================================
def load_summary_stats_only(base_dir="./outputs/summary"):

    all_data = []

    pattern = re.compile(
        r"summary_(?:([A-Za-z]+)_)?var(\d+)_deg([\d.]+)_samples(\d+)_alpha([\d.]+)_([A-Za-z0-9\-_]+)_variant(\d+)\.csv"
    )

    for root, _, files in os.walk(base_dir):
        for fname in files:

            if not fname.startswith("summary_") or not fname.endswith(".csv"):
                continue

            m = pattern.match(fname)
            if m is None:
                continue

            graph_type, var, deg, samples, alpha, algo, variant = m.groups()

            if graph_type is None:
                graph_type = "ER" if "ER" in root else "SF"

            if graph_type != "ER":
                continue

            file_path = os.path.join(root, fname)

            try:
                with open(file_path, "r") as f:
                    lines = f.readlines()

                summary_indices = [
                    i for i, line in enumerate(lines)
                    if "Summary (Mean/Std)" in line
                ]

                if not summary_indices:
                    continue

                start_idx = summary_indices[-1]
                summary_content = "".join(lines[start_idx + 1:])
                df_stats = pd.read_csv(StringIO(summary_content))

                if df_stats.columns[0].startswith("Unnamed"):
                    df_stats.rename(columns={df_stats.columns[0]: "metric"}, inplace=True)

                df_stats["metric"] = df_stats["metric"].astype(str).str.strip()
                df_stats = df_stats.drop_duplicates(subset=["metric"], keep="last")
                df_stats.set_index("metric", inplace=True)

                curr_dict = {
                    "GraphType": graph_type,
                    "Var": int(var),
                    "Deg": round(float(deg), 1),
                    "Samples": int(samples),
                    "Alpha": round(float(alpha), 1),
                    "Algorithm": algo.upper(),
                    "Variant": int(variant),
                }

                for metric in ["time", "ecshd", "benchpress_f1"]:
                    if metric in df_stats.index:
                        curr_dict[f"{metric}_mean"] = float(df_stats.loc[metric, "mean"])
                        curr_dict[f"{metric}_std"] = float(df_stats.loc[metric, "std"])
                    else:
                        curr_dict[f"{metric}_mean"] = np.nan
                        curr_dict[f"{metric}_std"] = np.nan

                all_data.append(curr_dict)

            except Exception as e:
                print(f"Error parsing {fname}: {e}")

    if not all_data:
        raise RuntimeError("❌ No summary stats found.")

    return pd.DataFrame(all_data)


# =========================================================
# 2. GENERATE SINGLE ALGORITHM SIDEWAYSTABLE
# =========================================================
def generate_single_algo_table(df, algo, d, output_dir):

    df = df[df["Algorithm"] == algo].copy()

    if df.empty:
        print(f"⚠ No data for {algo}")
        return

    target_variants = [0, 1, 2, 3]
    target_metrics = ["ecshd", "benchpress_f1", "time"]

    metric_dir = {
        "ecshd": "min",
        "benchpress_f1": "max",
        "time": "min",
    }

    grouped = df.groupby(["Var", "Deg", "Samples", "Alpha"])
    rows = []

    prev_deg = None

    for (var, deg, samples, alpha), group in grouped:

        if prev_deg is not None and deg != prev_deg:
            rows.append(r"\midrule")

        prev_deg = deg

        n_str = f"10^{int(np.log10(samples))}"
        row = f"{var} & {deg} & ${n_str}$ & {alpha}"


        rank_vals = {}

        for m in target_metrics:
            vals = []

            for v in target_variants:
                g = group[group["Variant"] == v]
                if not g.empty:
                    val = g[f"{m}_mean"].values[0]
                    if not pd.isna(val):
                        vals.append(val)

            if len(vals) == 0:
                rank_vals[m] = (None, None)
                continue

            vals_sorted = sorted(vals)
            if metric_dir[m] == "max":
                vals_sorted = vals_sorted[::-1]

            best = vals_sorted[0]
            second = vals_sorted[1] if len(vals_sorted) > 1 else None

            rank_vals[m] = (best, second)


        for m in target_metrics:

            vals = []

            if m == "time":
                ranking_variants = [v for v in target_variants if v != 0]
            else:
                ranking_variants = target_variants

            for v in ranking_variants:
                g = group[group["Variant"] == v]
                if not g.empty:
                    val = g[f"{m}_mean"].values[0]
                    if not pd.isna(val):
                        vals.append(val)

            if len(vals) == 0:
                rank_vals[m] = (None, None)
                continue

            vals_sorted = sorted(vals)
            if metric_dir[m] == "max":
                vals_sorted = vals_sorted[::-1]

            best = vals_sorted[0]
            second = vals_sorted[1] if len(vals_sorted) > 1 else None

            rank_vals[m] = (best, second)

            best_val, second_val = rank_vals[m]

            for v in target_variants:

                g = group[group["Variant"] == v]

                if g.empty:
                    row += " & -"
                    continue

                mean = g[f"{m}_mean"].values[0]
                std = g[f"{m}_std"].values[0]

                if pd.isna(mean):
                    row += " & -"
                    continue


                if m == "ecshd":
                    mean_str = f"{int(mean)}"
                    std_str = f"{std:.1f}"
                elif m == "benchpress_f1":
                    mean_str = f"{mean:.3f}"
                    std_str = f"{std:.2f}"
                else:
                    mean_str = f"{mean:.2f}"
                    std_str = f"{std:.2f}"

                pm_str = r"\boldsymbol{\pm}"
                base_cell = f"{mean_str}\\,{{\\scriptscriptstyle {pm_str}\\,{std_str}}}"
                # base_cell = f"{mean_str}{{({std_str})}}"

                if best_val is not None and abs(mean - best_val) < 1e-7:
                    cell = f"$\\boldsymbol{{{base_cell}}}$"

                elif second_val is not None and abs(mean - second_val) < 1e-7:
                    cell = f"$\\underline{{{base_cell}}}$"

                else:
                    cell = f"${base_cell}$"

                row += f" & {cell}"

        row += r" \\"   
        rows.append(row)

    # =====================================================
    # LATEX STRUCTURE
    # =====================================================
    latex = rf"""
\begin{{sidewaystable*}}[t]
\centering
\scriptsize
\renewcommand{{\arraystretch}}{{0.6}}
\setlength{{\tabcolsep}}{{2pt}}

\setlength{{\aboverulesep}}{{0.5pt}}
\setlength{{\belowrulesep}}{{0.5pt}}
\setlength{{\cmidrulesep}}{{0.5pt}}

\begin{{tabularx}}{{\textwidth}}{{@{{}}
c c c c |
>{{\raggedleft\arraybackslash}}X
>{{\raggedleft\arraybackslash}}X
>{{\columncolor{{gray!15}}\raggedleft\arraybackslash}}X
>{{\columncolor{{gray!15}}\raggedleft\arraybackslash}}X |
>{{\raggedleft\arraybackslash}}X
>{{\raggedleft\arraybackslash}}X
>{{\columncolor{{gray!15}}\raggedleft\arraybackslash}}X
>{{\columncolor{{gray!15}}\raggedleft\arraybackslash}}X |
>{{\raggedleft\arraybackslash}}X
>{{\raggedleft\arraybackslash}}X
>{{\columncolor{{gray!15}}\raggedleft\arraybackslash}}X
>{{\columncolor{{gray!15}}\raggedleft\arraybackslash}}X
@{{}}}}
\toprule

\multicolumn{{16}}{{c}}{{\bfseries \small {algo} Variants ($d={d}$)}} \\
\midrule

\multicolumn{{4}}{{c|}}{{Configuration}} &
\multicolumn{{4}}{{c|}}{{SHD ($\downarrow$)}} &
\multicolumn{{4}}{{c|}}{{F1 ($\uparrow$)}} &
\multicolumn{{4}}{{c}}{{RT ($\downarrow$)}} \\
\cmidrule(lr){{1-4}}
\cmidrule(lr){{5-8}}
\cmidrule(lr){{9-12}}
\cmidrule(lr){{13-16}}

\bfseries $d$ & \bfseries $\rho$ & \bfseries $n$ & \bfseries $\lambda$ &
\bfseries {algo} &
\bfseries {algo}-\textsc{{d}} &
\cellcolor{{gray!20}}\bfseries {algo}-\textsc{{dp}} &
\cellcolor{{gray!20}}\bfseries {algo}-\textsc{{dp+}} &
\bfseries {algo} &
\bfseries {algo}-\textsc{{d}} &
\cellcolor{{gray!20}}\bfseries {algo}-\textsc{{dp}} &
\cellcolor{{gray!20}}\bfseries {algo}-\textsc{{dp+}} &
\bfseries {algo} &
\bfseries {algo}-\textsc{{d}} &
\cellcolor{{gray!20}}\bfseries {algo}-\textsc{{dp}} &
\cellcolor{{gray!20}}\bfseries {algo}-\textsc{{dp+}} \\
\midrule
"""

    latex += "\n".join(rows)

    latex += rf"""
\bottomrule
\end{{tabularx}}
\label{{table:{algo.lower()}_grid}}
\end{{sidewaystable*}}
"""

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{algo.lower()}_grid.tex")

    with open(out_path, "w") as f:
        f.write(latex)

    print(f"✅ Saved: {out_path}") 


# =========================================================
# 3. RUN
# =========================================================
raw_stats = load_summary_stats_only("./outputs/summary")

filtered_data = raw_stats[
    raw_stats["Var"].isin([50, 100, 200]) &
    raw_stats["Deg"].isin([2.0, 3.0, 5.0]) &
    raw_stats["Samples"].isin([1000, 10000, 100000]) &
    raw_stats["Alpha"].isin([1.0, 2.0, 4.0])
]

# algos = sorted(filtered_data["Algorithm"].unique())

# for algo in algos:
#     generate_single_algo_table(
#         filtered_data,
#         algo,
#         "./outputs/tables/supplementary"
#     )

for d in [50, 100, 200]:

    df_d = filtered_data[filtered_data["Var"] == d]

    for algo in sorted(df_d["Algorithm"].unique()):
        generate_single_algo_table(
            df_d,
            algo,
            d,  
            f"./outputs/tables/supplementary/d{d}"
        )
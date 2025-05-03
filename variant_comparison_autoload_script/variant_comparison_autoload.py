import os
import re
import subprocess
import sys


# Absolute or relative path to the directory you want to include
env_dir = "/data/yaoscript"

# Add to sys.path if not already included
if env_dir not in sys.path:
    sys.path.append(env_dir)

from io import StringIO
from itertools import combinations

# Auto-install required packages if missing
def ensure_dependencies():
    required = ["pandas", "matplotlib", "seaborn", "scipy"]
    for pkg in required:
        try:
            __import__(pkg)
        except ImportError:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])




def load_csv_files(directory, group_name, report_type):
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)  # ← centralized path join
                if report_type == 'coverage' and 'CoverageRegion' in file:
                    file_paths.append(file_path)
                elif report_type == 'variant' and 'CoverageRegion' not in file:
                    file_paths.append(file_path)


    dataframes = []
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()

            header_line = None
            if report_type == 'coverage':
                header_line = next((l.strip().replace('\t', ',') for l in lines if l.startswith('#RegionName')), None)
            elif report_type == 'variant':
                header_line = next((l.strip().replace('\t', ',') for l in lines if l.startswith('#Gene')), None)

            if not header_line:
                continue

            header_cols = [col.strip() for col in header_line.split(',')]
            num_columns = len(header_cols)

            start_idx = next(i for i, line in enumerate(lines) if header_line in line.replace('\t', ','))
            valid_data_lines = [line for line in lines[start_idx:] if len(line.strip().split(',')) == num_columns]

            df = pd.read_csv(StringIO("".join(valid_data_lines)), sep=',', on_bad_lines='skip')

            if "HGVSCodingTranscript" in df.columns and "HGVSCoding" in df.columns:
                df["HGVSCoding"] = df["HGVSCodingTranscript"].astype(str) + ":" + df["HGVSCoding"].astype(str)

            if "#RegionName" in df.columns:
                df.rename(columns={"#RegionName": "RegionName"}, inplace=True)
            if "#Gene" in df.columns:
                df.rename(columns={"#Gene": "Gene"}, inplace=True)

            df.insert(0, "Group", group_name)
            sample_id = re.search(r'(MO\d{6})', os.path.basename(file_path))
            df.insert(0, "Sample ID", sample_id.group(1) if sample_id else os.path.basename(file_path))
            dataframes.append(df)

        except Exception as e:
            print(f"Error loading {file_path}: {e}")
            continue

    if not dataframes:
        raise ValueError(f"No valid CSV files found in {directory} for type {report_type}.")

    return pd.concat(dataframes, ignore_index=True)


# ==== MAIN SCRIPT ====
if __name__ == "__main__":
    ensure_dependencies()  # 🔥 Ensure this happens BEFORE plotting libraries are imported

    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from scipy.stats import linregress

    script_dir = script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    group_dirs = [
        os.path.join(script_dir, d)
        for d in os.listdir(script_dir)
        if os.path.isdir(os.path.join(script_dir, d)) and not d.startswith('.') and d != "output"
    ]
    group_names = [os.path.basename(d) for d in group_dirs]
    group_info = list(zip(group_names, group_dirs))

    # Load data per group
    group_cov_dfs = {
        group_name: load_csv_files(path, group_name, report_type='coverage')
        .assign(AverageCoverage=lambda df: pd.to_numeric(df["AverageCoverage"], errors='coerce')
        )
        for group_name, path in group_info
    }

    group_var_dfs = {
        group_name: load_csv_files(path, group_name, report_type='variant')
        for group_name, path in group_info
    }

    # Coverage analysis (line plot per group)
    available_cov = [g for g in group_cov_dfs if 'AverageCoverage' in group_cov_dfs[g].columns and 'RegionName' in group_cov_dfs[g].columns]
    if available_cov:
        coverage_df = pd.concat([
            group_cov_dfs[g][["RegionName", "AverageCoverage"]].assign(Group=g)
            for g in available_cov
        ])
        coverage_df = coverage_df.groupby(["RegionName", "Group"]).mean().reset_index()
        pivot_cov = coverage_df.pivot(index="RegionName", columns="Group", values="AverageCoverage").fillna(0)
        pivot_cov = pivot_cov.sort_index()

        plt.figure(figsize=(20, 6))
        for group in pivot_cov.columns:
            plt.plot(pivot_cov.index, pivot_cov[group], marker='o', label=group)

        plt.title("Average Coverage per RegionName across Groups")
        plt.ylabel("Average Coverage")
        plt.xlabel("RegionName")
        plt.legend()
        plt.tight_layout()
        plt.xticks([])
        ax = plt.gca()
        for i, label in enumerate(pivot_cov.index):
            # X-axis position for the tick/label
            x_pos = i
            # Alternate label vertical position
            y_offset = -0.15 if i % 2 == 0 else -0.09
            # Draw connecting vertical line from axis to label
            ax.plot([x_pos, x_pos], [0, y_offset + 0.01], transform=ax.get_xaxis_transform(),
                    color='gray', lw=0.5, linestyle='-')
            # Draw label
            ax.text(x_pos, y_offset, label, ha='right', va='top', fontsize=5, rotation=90,
                    transform=ax.get_xaxis_transform())

        # Limit x range and adjust layout
        plt.xlim(-1, len( pivot_cov .index))
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.3)
        plt.savefig(os.path.join(output_dir, "coverage_comparison.png"))
        plt.close()


        print("✅ Coverage comparison plot saved as 'output/coverage_comparison.png'")

    # Variant concordance analysis (pairwise comparison)
    # Variant concordance analysis (pairwise comparison)
    merge_keys = ["Sample ID", "Gene", "HGVSCoding","Chr:ChrPos"]
    metrics = ["VariantFrequency", "Coverage", "ReadBalance", "Pathogenicity", "Type"]

    comparison_rows = []
    for group1, group2 in combinations(group_var_dfs.keys(), 2):
        df1 = group_var_dfs[group1]
        df2 = group_var_dfs[group2]

        if not all(k in df1.columns for k in merge_keys) or not all(k in df2.columns for k in merge_keys):
            continue

        df1 = df1[merge_keys + metrics].copy()
        df2 = df2[merge_keys + metrics].copy()

        df_merged = df1.merge(df2, on=merge_keys, how='outer', suffixes=(f'_{group1}', f'_{group2}'))
        df_merged.fillna(0, inplace=True)

        for _, row in df_merged.iterrows():
            v1 = row.get(f"VariantFrequency_{group1}", 0)
            v2 = row.get(f"VariantFrequency_{group2}", 0)
            percent_diff = abs(v1 - v2) / max(v1, v2) if max(v1, v2) != 0 else 0

            record = {
                "Sample ID": row["Sample ID"],
                "Gene": row["Gene"],
                "HGVSCoding": row["HGVSCoding"],
                "Chr:ChrPos":row["Chr:ChrPos"],
                "Group_1": group1,
                "Group_2": group2,
                "VariantFrequency_diff_gt_50%": percent_diff > 0.5,
            }

            for metric in metrics:
                record[f"{metric}_{group1}"] = row.get(f"{metric}_{group1}", None)
                record[f"{metric}_{group2}"] = row.get(f"{metric}_{group2}", None)

            comparison_rows.append(record)

    comparison_df = pd.DataFrame(comparison_rows)

    # Sort the columns to group metrics of each group together
    id_cols = ["Sample ID", "Gene", "HGVSCoding", "Chr:ChrPos","Group_1", "Group_2", "VariantFrequency_diff_gt_50%"]
    metric_cols = []
    for metric in metrics:
        for group1, group2 in combinations(group_var_dfs.keys(), 2):
            g1 = f"{metric}_{group1}"
            g2 = f"{metric}_{group2}"
            if g1 not in metric_cols:
                metric_cols.append(g1)
            if g2 not in metric_cols:
                metric_cols.append(g2)

    final_cols = id_cols + [col for col in metric_cols if col in comparison_df.columns]
    comparison_df = comparison_df[final_cols]

    comparison_df.to_csv(os.path.join(output_dir, "variant_group_pairwise_comparison.csv"), index=False)
    print("✅ Exported variant comparison table to 'output/variant_group_pairwise_comparison.csv'")

    # Correlation plots
    for group1, group2 in combinations(group_var_dfs.keys(), 2):
        df1 = group_var_dfs[group1][merge_keys + ["VariantFrequency"]].copy()
        df2 = group_var_dfs[group2][merge_keys + ["VariantFrequency"]].copy()
        df_corr = df1.merge(df2, on=merge_keys, how='outer', suffixes=(f'_{group1}', f'_{group2}')).fillna(0)
        # Drop rows where both VariantFrequency values are NaN
        df_corr = df_corr.dropna(subset=[f"VariantFrequency_{group1}", f"VariantFrequency_{group2}"], how='all')

        x = df_corr[f"VariantFrequency_{group1}"].fillna(0)
        y = df_corr[f"VariantFrequency_{group2}"].fillna(0)
        plt.figure(figsize=(6, 6))
        sns.scatterplot(x=x, y=y)

        # Fit regression line
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = slope * x + intercept
        plt.plot(x, line, color='red', linestyle='--')
        # Annotate R² and regression equation
        plt.text(0.05, 0.9, f'R² = {r_value ** 2:.3f}', transform=plt.gca().transAxes, fontsize=10)
        plt.text(0.05, 0.85, f'y = {slope:.2f}x + {intercept:.2f}', transform=plt.gca().transAxes, fontsize=10)
        plt.title(f'{group1} vs {group2} Correlation')
        plt.grid(True)
        plt.tight_layout()
        plt.xlabel(f"{group1} VariantFrequency")
        plt.ylabel(f"{group2} VariantFrequency")
        plt.savefig(os.path.join(output_dir, f"correlation_{group1}_vs_{group2}.png"))
        plt.close()

    print("✅ Correlation plots saved for all group combinations")

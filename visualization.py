import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from config import BASE_WORK_DIR, MUTATION_RATES

def collect_all_results(base_dir):
    data = []
    for mut_rate in MUTATION_RATES:
        work_dir = f"{base_dir}_{mut_rate}"
        if not os.path.isdir(work_dir):
            continue
        for rep_dir in os.listdir(work_dir):
            rep_path = os.path.join(work_dir, rep_dir)
            if not os.path.isdir(rep_path):
                continue
            if not (rep_dir.startswith("m") and "_sel" in rep_dir and "_rep" in rep_dir):
                continue
            result_file = os.path.join(rep_path, "result.txt")
            if not os.path.isfile(result_file):
                continue
            parts = rep_dir.split("_")
            try:
                m = float(parts[0][1:])
                n_sites = int(parts[1][3:])
                rep = int(parts[2][3:])
            except Exception as e:
                continue
            with open(result_file, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    rf_line = lines[0].strip()
                    wrf_line = lines[1].strip()
                    if rf_line.startswith("RF="):
                        rf = float(rf_line[3:])
                    else:
                        rf = float(rf_line)
                    if wrf_line.startswith("WRF="):
                        wrf = float(wrf_line[4:])
                    else:
                        wrf = float(wrf_line)
                    data.append({
                        "mutation_rate": mut_rate,
                        "m": m,
                        "n_selected_sites": n_sites,
                        "replicate": rep,
                        "RF": rf,
                        "WRF": wrf
                    })
    return pd.DataFrame(data)

def create_tables_and_plots(df, metric="RF"):
    if df.empty:
        return
    summary = df.groupby(['mutation_rate', 'm', 'n_selected_sites'])[metric].agg(['mean', 'std', 'min', 'max']).round(6)
    summary.columns = ['mean', 'std', 'min', 'max']
    summary = summary.reset_index()
    for mut_rate in df['mutation_rate'].unique():
        subset = df[df['mutation_rate'] == mut_rate]
        pivot = subset.pivot_table(index='n_selected_sites', columns='m', values=metric, aggfunc='mean')
        if pivot.empty:
            continue
        plt.figure(figsize=(8,6))
        sns.heatmap(pivot, annot=True, fmt='.4f', cmap='YlOrRd', cbar_kws={'label': metric})
        plt.title(f'{metric} heatmap (mutation_rate={mut_rate})')
        plt.xlabel('m')
        plt.ylabel('n_selected_sites')
        plt.tight_layout()
        plt.savefig(f'heatmap_{metric}_mut{mut_rate}.png', dpi=150)
        plt.close()
    plt.figure(figsize=(10,6))
    for mut_rate in df['mutation_rate'].unique():
        for n_sites in sorted(df['n_selected_sites'].unique()):
            sub = df[(df['mutation_rate']==mut_rate) & (df['n_selected_sites']==n_sites)]
            means = sub.groupby('m')[metric].mean()
            if means.empty:
                continue
            plt.plot(means.index, means.values, marker='o', label=f'mut={mut_rate}, sites={n_sites}')
    plt.xlabel('m')
    plt.ylabel(metric)
    plt.title(f'{metric} vs m')
    plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'lineplot_{metric}_vs_m.png', dpi=150)
    plt.close()

def main():
    df = collect_all_results(BASE_WORK_DIR)
    if df.empty:
        return
    df.to_csv("all_results.csv", index=False)
    create_tables_and_plots(df, metric="RF")
    create_tables_and_plots(df, metric="WRF")
if __name__ == "__main__":
    main()
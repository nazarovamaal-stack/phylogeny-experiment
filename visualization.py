import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from config import BASE_WORK_DIR, MUTATION_RATES

def collect_all_results(base_dir):
    """Собирает все файлы result.txt из папок exp_stable_*/m*_sel*_rep*/ и возвращает DataFrame с колонками mutation_rate, m, n_sites, replicate, RF, WRF."""
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

def mann_whitney_test(df, results_dir):
    """Выполняет тест Манна-Уитни (контроль m=0 vs m>0)."""
    metrics = ['RF', 'WRF']
    mut_rates = sorted(df['mutation_rate'].unique())
    n_sites_vals = sorted(df['n_selected_sites'].unique())
    results = []
    for mut_rate in mut_rates:
        for n_sites in n_sites_vals:
            subset = df[(df['mutation_rate'] == mut_rate) & (df['n_selected_sites'] == n_sites)]
            if subset.empty or 0 not in subset['m'].values:
                continue
            control = subset[subset['m'] == 0]
            for metric in metrics:
                control_vals = control[metric].values
                for m_val in sorted([m for m in subset['m'].unique() if m > 0]):
                    test_vals = subset[subset['m'] == m_val][metric].values
                    if len(control_vals) == 0 or len(test_vals) == 0:
                        continue
                    stat, p = mannwhitneyu(control_vals, test_vals, alternative='less')
                    results.append({
                        'mutation_rate': mut_rate,
                        'n_selected_sites': n_sites,
                        'metric': metric,
                        'm_control': 0,
                        'm_test': m_val,
                        'U_statistic': stat,
                        'p_value': p
                    })
    df_results = pd.DataFrame(results)
    df_results.to_csv(os.path.join(results_dir, 'mann_whitney_all.csv'), index=False)

def create_tables_and_plots(df, results_dir, metric="RF"):
    """По данным df строит тепловые карты и линейные графики зависимости метрики от m."""
    if df.empty:
        return
    heatmap_dir = os.path.join(results_dir, 'heatmaps')
    os.makedirs(heatmap_dir, exist_ok=True)
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
        plt.savefig(os.path.join(heatmap_dir, f'heatmap_{metric}_mut{mut_rate}.png'), dpi=150)
        plt.close()
    lineplot_dir = os.path.join(results_dir, 'lineplots')
    os.makedirs(lineplot_dir, exist_ok=True)
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
    plt.savefig(os.path.join(lineplot_dir, f'lineplot_{metric}_vs_m.png'), dpi=150)
    plt.close()

def create_boxplots(df, results_dir):
    """Строит boxplot-графики для RF и WRF по m для каждой комбинации mutation_rate и n_selected_sites."""
    boxplot_dir = os.path.join(results_dir, 'boxplots')
    os.makedirs(boxplot_dir, exist_ok=True)
    for metric in ['RF', 'WRF']:
        for mut_rate in df['mutation_rate'].unique():
            for n_sites in sorted(df['n_selected_sites'].unique()):
                subset = df[(df['mutation_rate'] == mut_rate) & (df['n_selected_sites'] == n_sites)]
                if subset.empty:
                    continue
                plt.figure(figsize=(8,6))
                sns.boxplot(data=subset, x='m', y=metric)
                plt.title(f'{metric} | mutation_rate={mut_rate}, n_sites={n_sites}')
                plt.xlabel('m')
                plt.ylabel(metric)
                plt.grid(axis='y', alpha=0.3)
                plt.tight_layout()
                plt.savefig(os.path.join(boxplot_dir, f'boxplot_{metric}_mut{mut_rate}_n{n_sites}.png'), dpi=150)
                plt.close()

def main():
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    df = collect_all_results(BASE_WORK_DIR)
    if df.empty:
        return
    df.to_csv(os.path.join(results_dir, 'all_results.csv'), index=False)
    create_tables_and_plots(df, results_dir, metric="RF")
    create_tables_and_plots(df, results_dir, metric="WRF")
    create_boxplots(df, results_dir)
    mann_whitney_test(df, results_dir)

if __name__ == "__main__":
    main()
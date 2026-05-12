"""Визуализация результатов эксперимента"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def create_results_tables(results, output_dir, metric="metric"):
    """Создаёт сводные таблицы для заданной метрики."""
    df = pd.DataFrame(results)
    if df.empty:
        return df, None

    summary = df.groupby(['m', 'n_selected_sites'])['distance'].agg(['mean', 'std', 'min', 'max']).round(4)
    summary.columns = ['mean', 'std', 'min', 'max']
    summary = summary.reset_index()

    tables_dir = os.path.join(output_dir, f"tables_{metric}")
    os.makedirs(tables_dir, exist_ok=True)

    df.to_csv(os.path.join(tables_dir, "raw_results.csv"), index=False)
    summary.to_csv(os.path.join(tables_dir, "summary_results.csv"), index=False)

    # Markdown
    md = f"# Результаты ({metric})\n\n## Сводная статистика\n\n"
    md += "| m | Сайты | mean | std | min | max |\n|----|-------|------|-----|-----|-----|\n"
    for _, row in summary.iterrows():
        md += f"| {row['m']:.2f} | {int(row['n_selected_sites'])} | {row['mean']:.4f} | {row['std']:.4f} | {row['min']:.0f} | {row['max']:.0f} |\n"
    with open(os.path.join(tables_dir, "results_table.md"), 'w', encoding='utf-8') as f:
        f.write(md)

    # LaTeX
    latex = "\\begin{table}[h]\n\\centering\n\\caption{Сводные результаты (" + metric + ")}\n"
    latex += "\\begin{tabular}{ccccc}\n\\hline\nm & Сайты & Среднее & Станд.откл. & Диапазон \\\\\n\\hline\n"
    for _, row in summary.iterrows():
        latex += f"{row['m']:.2f} & {int(row['n_selected_sites'])} & {row['mean']:.4f} & {row['std']:.4f} & [{row['min']:.0f}, {row['max']:.0f}] \\\\\n"
    latex += "\\hline\n\\end{tabular}\n\\end{table}\n"
    with open(os.path.join(tables_dir, "results_table.tex"), 'w', encoding='utf-8') as f:
        f.write(latex)

    print(f"Таблицы для {metric} сохранены в {tables_dir}")
    return df, summary
def create_all_plots(results, output_dir, metric="metric"):
    """Создаёт все графики для одной метрики."""
    df = pd.DataFrame(results)
    if df.empty:
        return
    plots_dir = os.path.join(output_dir, f"plots_{metric}")
    os.makedirs(plots_dir, exist_ok=True)
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    create_heatmap(df, plots_dir, metric)
    create_line_plots(df, plots_dir, metric)
    create_box_plot(df, plots_dir, metric)
    create_distribution_plot(df, plots_dir, metric)
    create_3d_surface(df, plots_dir, metric)
    create_violin_plot(df, plots_dir, metric)

    print(f"Графики для {metric} сохранены в {plots_dir}")


def create_heatmap(df, output_dir, metric):
    pivot = df.pivot_table(values='distance', index='n_selected_sites', columns='m', aggfunc='mean')
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': f'{metric} расстояние'}, ax=ax)
    ax.set_title(f'Тепловая карта средних {metric}-расстояний', fontsize=14, fontweight='bold')
    ax.set_xlabel('Коэффициент отбора (m)', fontsize=12)
    ax.set_ylabel('Количество сайтов под отбором', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'heatmap_{metric}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def create_line_plots(df, output_dir, metric):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for n_sel in df['n_selected_sites'].unique():
        data = df[df['n_selected_sites'] == n_sel].groupby('m')['distance'].agg(['mean', 'std'])
        axes[0].errorbar(data.index, data['mean'], yerr=data['std'], marker='o', label=f'{n_sel} сайтов', capsize=5)
    axes[0].set_xlabel('Коэффициент отбора (m)')
    axes[0].set_ylabel(f'{metric} расстояние')
    axes[0].set_title(f'Зависимость {metric} от m')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    for m in df['m'].unique():
        data = df[df['m'] == m].groupby('n_selected_sites')['distance'].agg(['mean', 'std'])
        axes[1].errorbar(data.index, data['mean'], yerr=data['std'], marker='s', label=f'm={m}', capsize=5)
    axes[1].set_xlabel('Количество сайтов под отбором')
    axes[1].set_ylabel(f'{metric} расстояние')
    axes[1].set_title(f'Зависимость {metric} от числа сайтов')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'line_plots_{metric}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def create_box_plot(df, output_dir, metric):
    fig, ax = plt.subplots(figsize=(12, 6))
    df['condition'] = df.apply(lambda x: f"m={x['m']}, s={x['n_selected_sites']}", axis=1)
    df.boxplot(column='distance', by='condition', ax=ax, rot=45)
    ax.set_title(f'Распределение {metric}-расстояний по условиям', fontsize=14, fontweight='bold')
    ax.set_xlabel('Условия эксперимента')
    ax.set_ylabel(f'{metric} расстояние')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'box_plot_{metric}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def create_distribution_plot(df, output_dir, metric):
    fig, ax = plt.subplots(figsize=(10, 6))
    for m in df['m'].unique():
        data = df[df['m'] == m]['distance']
        ax.hist(data, alpha=0.5, label=f'm={m}', bins=15)
    ax.set_xlabel(f'{metric} расстояние')
    ax.set_ylabel('Частота')
    ax.set_title(f'Распределение {metric}-расстояний')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'distribution_{metric}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def create_3d_surface(df, output_dir, metric):
    from mpl_toolkits.mplot3d import Axes3D
    pivot = df.pivot_table(values='distance', index='n_selected_sites', columns='m', aggfunc='mean')
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(pivot.columns, pivot.index)
    Z = pivot.values
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8, linewidth=0.5, edgecolor='black')
    ax.set_xlabel('Коэффициент отбора (m)')
    ax.set_ylabel('Сайты под отбором')
    ax.set_zlabel(f'{metric} расстояние')
    ax.set_title(f'3D поверхность {metric}-расстояний')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'3d_surface_{metric}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def create_violin_plot(df, output_dir, metric):
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(data=df, x='n_selected_sites', y='distance', hue='m', split=True, ax=ax, palette='Set2')
    ax.set_xlabel('Количество сайтов под отбором')
    ax.set_ylabel(f'{metric} расстояние')
    ax.set_title(f'Violin plot распределения {metric}-расстояний')
    ax.legend(title='m')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'violin_{metric}.png'), dpi=300, bbox_inches='tight')
    plt.close()
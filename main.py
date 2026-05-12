"""Основной скрипт эксперимента"""

import os
import csv
import shutil
from collections import defaultdict
from config import (
    M_VALUES, NUM_SELECTED_SITES, N_TAXA,
    NEUTRAL_LENGTH, N_REPLICATES, WORK_DIR
)
from utils.system import ensure_dir
from utils.simulation import run_vgsim_wrapper, run_alisim
from utils.sequence import concatenate_fasta
from utils.phylogeny import run_iqtree, rf_distance, branch_score_distance
from utils.visualization import create_results_tables, create_all_plots


def main():
    """Главная функция эксперимента"""
    ensure_dir(WORK_DIR)
    results_rf = []
    results_bsd = []
    n_comb = len(M_VALUES) * len(NUM_SELECTED_SITES) * N_REPLICATES
    print(f"количество комбинаций: {n_comb}\n")
    for m in M_VALUES:
        for n_sel in NUM_SELECTED_SITES:
            rf_values = []
            bsd_values = []
            for rep in range(1, N_REPLICATES + 1):
                rep_dir = os.path.join(WORK_DIR, f"m{m}_sel{n_sel}_rep{rep}")
                ensure_dir(rep_dir)
                seed = 42 + rep * 100 + n_sel * 1000 + int(m * 100)
                print(f"m={m}, sites={n_sel}, rep={rep}")
                try:
                    if n_sel > 0:
                        selected_fasta = os.path.join(rep_dir, "selected.fasta")
                        true_tree_path = run_vgsim_wrapper(
                            N_TAXA, n_sel, m, selected_fasta, seed
                        )
                    else:
                        selected_fasta = os.path.join(rep_dir, "selected.fasta")
                        true_tree_path = run_vgsim_wrapper(
                            N_TAXA, 0, m, selected_fasta, seed
                        )
                        if os.path.exists(selected_fasta):
                            os.remove(selected_fasta)
                    neutral_fasta = os.path.join(rep_dir, "neutral.fasta")
                    run_alisim(true_tree_path, neutral_fasta, seed)
                    combined_fasta = os.path.join(rep_dir, "combined.fasta")
                    if n_sel > 0:
                        concatenate_fasta(neutral_fasta, selected_fasta, combined_fasta)
                    else:
                        shutil.copy(neutral_fasta, combined_fasta)
                    iqtree_prefix = os.path.join(rep_dir, "iqtree_out")
                    inferred_tree = run_iqtree(combined_fasta, iqtree_prefix, seed)
                    rf = rf_distance(true_tree_path, inferred_tree)
                    bsd = branch_score_distance(true_tree_path, inferred_tree)

                    results_rf.append({
                        "m": m,
                        "n_selected_sites": n_sel,
                        "replicate": rep,
                        "distance": rf
                    })
                    results_bsd.append({
                        "m": m,
                        "n_selected_sites": n_sel,
                        "replicate": rep,
                        "distance": bsd
                    })
                    rf_values.append(rf)
                    bsd_values.append(bsd)
                    print(f"RF={rf:.4f}")
                    print(f"BSD={bsd:.4f}")
                except Exception as e:
                    print(f"Ошибка в репликации: {e}")
                    continue
            if rf_values:
                print(f"\nИтоговое среднее для m={m}, sites={n_sel}:")
                print(f"   RF = {sum(rf_values)/len(rf_values):.4f}")
                print(f"   BSD = {sum(bsd_values)/len(bsd_values):.4f}")
                print('-' * 50)
    for metric, results in [("rf", results_rf), ("bsd", results_bsd)]:
        csv_path = os.path.join(WORK_DIR, f"results_{metric}.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["m", "n_selected_sites", "replicate", "distance"]
            )
            writer.writeheader()
            writer.writerows(results)
    create_results_tables(results_rf, WORK_DIR, metric="RF")
    create_results_tables(results_bsd, WORK_DIR, metric="BSD")
    create_all_plots(results_rf, WORK_DIR, metric="RF")
    create_all_plots(results_bsd, WORK_DIR, metric="BSD")
    print(f"Результаты в папке {WORK_DIR}")
if __name__ == "__main__":
    main()
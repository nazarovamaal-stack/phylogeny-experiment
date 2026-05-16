import os
import sys
from tqdm import tqdm
from config import M_VALUES, NUM_SELECTED_SITES, N_TAXA, NEUTRAL_LENGTH, N_REPLICATES, BASE_WORK_DIR, MUTATION_RATES
from utils.system import ensure_dir
from utils.simulation import run_vgsim_wrapper, run_alisim
from utils.phylogeny import run_iqtree, rf_distance, weighted_rf_distance

def save_result(rep_dir, m, n_sel, rep, rf, wrf):
    """Сохраняет RF и WRF в файл вида mX_selY_repZ.txt"""
    result_file = os.path.join(rep_dir, f"result.txt")
    with open(result_file, 'w') as f:
        f.write(f"RF={rf:.6f}\n")
        f.write(f"WRF={wrf:.6f}\n")

def run_experiment(mutation_rate, work_dir):
    """Запускает полный цикл эксперимента для заданной скорости мутаций: по всем m, n_sites, репликам."""
    ensure_dir(work_dir)
    n_comb = len(M_VALUES) * len(NUM_SELECTED_SITES) * N_REPLICATES
    print(f"\nMUTATION RATE = {mutation_rate}")
    print(f"Total combinations: {n_comb}\n")
    sys.stdout.flush()
    for m in tqdm(M_VALUES):
        for n_sel in NUM_SELECTED_SITES:
            for rep in range(1, N_REPLICATES + 1):
                rep_dir = os.path.join(work_dir, f"m{m}_sel{n_sel}_rep{rep}")
                result_file = os.path.join(rep_dir, "result.txt")
                if os.path.exists(result_file):
                    continue
                ensure_dir(rep_dir)
                seed = 42 + rep * 100 + n_sel * 1000 + int(m * 100)
                print(f"m={m}, sites={n_sel}, rep={rep}")
                sys.stdout.flush()
                effective_sites = n_sel if m > 0 else 0
                true_tree_path = None
                try:
                    dummy_fasta = os.path.join(rep_dir, "dummy.fasta")
                    true_tree_path = run_vgsim_wrapper(N_TAXA, effective_sites, m, dummy_fasta, seed, mutation_rate)
                except Exception as e:
                    sys.stdout.flush()
                if true_tree_path is None:
                    sys.stdout.flush()
                    continue
                try:
                    neutral_fasta = os.path.join(rep_dir, "neutral.fasta")
                    _, scaled_tree_path = run_alisim(true_tree_path, neutral_fasta, seed, mutation_rate, seq_len=NEUTRAL_LENGTH)
                    iqtree_prefix = os.path.join(rep_dir, "iqtree_out")
                    inferred_tree = run_iqtree(neutral_fasta, iqtree_prefix, seed)
                    rf = rf_distance(scaled_tree_path, inferred_tree)
                    wrf = weighted_rf_distance(scaled_tree_path, inferred_tree)
                    save_result(rep_dir, m, n_sel, rep, rf, wrf)
                    print(f"RF={rf:.6f}  WRF={wrf:.6f}")
                    sys.stdout.flush()
                except Exception as e:
                    sys.stdout.flush()
                    continue

def main():
    for mut_rate in MUTATION_RATES:
        work_dir = f"{BASE_WORK_DIR}_{mut_rate}"
        run_experiment(mut_rate, work_dir)

if __name__ == "__main__":
    main()
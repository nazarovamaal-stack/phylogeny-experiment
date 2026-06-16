import os
import sys
import numpy as np
import dendropy
import VGsim
import subprocess
import itertools
from config import (
    TRANSMISSION_RATE, RECOVERY_RATE, SAMPLING_RATE,
    SUBSTITUTION_WEIGHTS, POPULATION_SIZE,
    EPIDEMIC_TIME, MAX_ITERATIONS, NEUTRAL_LENGTH
)

def run_vgsim_wrapper(n_taxa, n_sites, m, out_fasta, seed, mutation_rate):
    """Генерирует истинное филогенетическое дерево с помощью VGsim с возможностью отбора."""
    actual_sites = max(1, n_sites)
    max_attempts = 5
    for attempt in range(max_attempts):
        current_seed = seed + attempt
        out_base = os.path.join(os.path.dirname(out_fasta), f"vgsim_tree_attempt{attempt}")
        try:
            simulator = VGsim.Simulator(actual_sites, 1, 1, seed=current_seed)
            simulator.set_transmission_rate(TRANSMISSION_RATE)
            simulator.set_recovery_rate(RECOVERY_RATE)
            simulator.set_sampling_rate(SAMPLING_RATE)
            simulator.set_mutation_rate(mutation_rate, SUBSTITUTION_WEIGHTS)
            if n_sites > 0 and m > 0:
                nucleotides = ['A', 'T', 'C', 'G']
                for haplotype in itertools.product(nucleotides, repeat=n_sites):
                    haplotype_str = ''.join(haplotype)
                    g_count = haplotype_str.count('G')
                    if g_count == 0:
                        rate = TRANSMISSION_RATE
                    else:
                        rate = TRANSMISSION_RATE * ((1.0 + m) ** g_count)
                    simulator.set_transmission_rate(rate, haplotype=haplotype_str)
            simulator.set_population_size(POPULATION_SIZE, population=0)
            simulator.simulate(MAX_ITERATIONS, n_taxa, EPIDEMIC_TIME, 'tau')
            simulator.genealogy()
            simulator.export_newick(out_base)
            tree_path = out_base + "_tree.nwk"
            if not os.path.exists(tree_path):
                alt_path = out_base + ".nwk"
                if os.path.exists(alt_path):
                    tree_path = alt_path
                else:
                    continue
            tree = dendropy.Tree.get_from_path(tree_path, schema="newick", rooting="force-unrooted")
        except Exception as e:
            sys.stdout.flush()
            continue
        leaves = list(tree.leaf_node_iter())
        if len(leaves) >= n_taxa:
            for i, leaf in enumerate(leaves[:n_taxa], start=1):
                leaf.taxon.label = f"taxon_{i}"
            if len(leaves) > n_taxa:
                tree.retain_taxa_with_labels([f"taxon_{i}" for i in range(1, n_taxa + 1)])
            tree.write(path=tree_path, schema="newick", suppress_rooting=True, unquoted_underscores=True)
            return tree_path
        else:
            sys.stdout.flush()
    sys.stdout.flush()
    return None


def run_alisim(scaled_tree_path, seq_len, out_base):
    """
    Запускает симулятор AliSim через командную строку IQ-TREE. Генерирует нейтральные последовательности по модели JC69.
    """
    out_fasta = out_base + ".fa"
    cmd = [
        "iqtree2",
        "--alisim", out_base,
        "-m", "JC",
        "-t", scaled_tree_path,
        "--num-sites", str(seq_len),
        "-af", "fasta",
        "--overwrite"
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.exists(out_fasta):
            return out_fasta
        else:
            raise FileNotFoundError(f"AliSim не создал файл: {out_fasta}")

    except subprocess.CalledProcessError as e:
        print(f"Ошибка при запуске AliSim: {e}")
        return None
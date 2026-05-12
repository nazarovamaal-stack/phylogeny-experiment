"""Утилиты для симуляции данных"""

import os
import sys
import numpy as np
import dendropy
from .system import run_cmd
from dendropy.model.discrete import Jc69, DiscreteCharacterEvolver
from config import IQTREE_CMD, VGSIM_CMD, NEUTRAL_LENGTH


def run_vgsim_wrapper(n_taxa, n_sites, m, out_fasta, seed):
    """
    Вызывает VGsim, возвращает путь к сгенерированному дереву (true_tree).
    Гарантирует, что дерево будет содержать ровно n_taxa таксонов.
    """
    if not os.path.exists(VGSIM_CMD):
        raise FileNotFoundError(f"VGsim_cmd.py не найден: {VGSIM_CMD}")

    rep_dir = os.path.dirname(os.path.abspath(out_fasta))
    base_name = os.path.join(rep_dir, "vgsim_run")
    tree_out = base_name + "_tree.nwk"
    mut_out = base_name + "_mutations.tsv"
    max_attempts = 100
    current_seed = seed

    for attempt in range(1, max_attempts + 1):
        cmd = [
            sys.executable, VGSIM_CMD,
            "--sampleSize", str(n_taxa),
            "--iterations", "1000",
            "--time", "-1",
            "--seed", str(current_seed),
            "--createNewick", base_name,
            "--writeMutations", base_name
        ]
        res = run_cmd(cmd)
        if not os.path.exists(tree_out):
            print("ошибка")

        tree = dendropy.Tree.get_from_path(
            tree_out, schema="newick", rooting="force-unrooted"
        )
        actual_taxa = len(list(tree.leaf_node_iter()))


        if actual_taxa == n_taxa:
            print(f"VGsim создал дерево с {actual_taxa} таксонами")
            break
        else:
            current_seed += 1
            if os.path.exists(tree_out):
                os.remove(tree_out)
            if os.path.exists(mut_out):
                os.remove(mut_out)
    else:
        print(f"Не удалось получить дерево с {n_taxa} таксонами за {max_attempts} попыток.")
        sys.exit(1)
    for i, leaf in enumerate(tree.leaf_node_iter(), 1):
        leaf.taxon.label = f"taxon_{i}"

    SCALE_FACTOR = 2.0
    for node in tree.preorder_node_iter():
        if node.edge_length is not None:
            node.edge_length /= SCALE_FACTOR
    tree.write(path=tree_out, schema="newick",
               suppress_rooting=True, unquoted_underscores=True)
    if n_sites > 0:
        _apply_selection_model(tree, n_sites, m, current_seed, out_fasta)
    else:
        with open(out_fasta, "w") as f:
            for leaf in tree.leaf_node_iter():
                f.write(f">{leaf.taxon.label}\n\n")

    if os.path.exists(mut_out):
        os.remove(mut_out)

    return tree_out

def _apply_selection_model(tree, n_sites, m, seed, out_fasta):
    """
    4-буквенная конвергенция: две фиксированные группы по 15 листьев.
    Для каждого сайта случайно выбирается пара разных нуклеотидов.
    При вероятности m/(1+m) одна группа получает первый нуклеотид, другая – второй,
    иначе обе группы получают одинаковый случайный нуклеотид (нейтрально).
    """
    import random
    rng = random.Random(seed)

    all_leaves = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    seqs = {name: [] for name in all_leaves}
    shuffled = all_leaves[:]
    rng.shuffle(shuffled)
    mid = len(shuffled) // 2
    group1 = set(shuffled[:mid])
    group2 = set(shuffled[mid:])
    prob_conv = m / (1.0 + m)
    bases = ['A', 'C', 'G', 'T']
    for _ in range(n_sites):
        if rng.random() < prob_conv:
            b1, b2 = rng.sample(bases, 2)
            for name in all_leaves:
                seqs[name].append(b1 if name in group1 else b2)
        else:
            b = rng.choice(bases)
            for name in all_leaves:
                seqs[name].append(b)

    with open(out_fasta, 'w') as f:
        for name in all_leaves:
            f.write(f">{name}\n{''.join(seqs[name])}\n")

def run_alisim(tree_path, out_fasta, seed, seq_len=None):
    """
    Нейтральная симуляция на том же дереве (модель JC69).
    Использует прямой обход дерева, без внешних программ.
    Гарантирует, что последовательности эволюционируют строго на переданном дереве.
    """
    if seq_len is None:
        seq_len = NEUTRAL_LENGTH
    import numpy as np
    tree = dendropy.Tree.get_from_path(
        tree_path, schema="newick", rooting="force-unrooted"
    )
    rng = np.random.default_rng(seed)
    node_to_seq = {}

    def simulate_node(node, parent_seq=None):
        if parent_seq is None:
            seq = rng.choice(['A', 'C', 'G', 'T'], size=seq_len).tolist()
        else:
            t = node.edge_length if node.edge_length is not None else 0.0
            if t > 0:
                p_same = 0.25 + 0.75 * np.exp(-4.0/3.0 * t)
                p_diff = 0.25 - 0.25 * np.exp(-4.0/3.0 * t)
                parent_arr = np.array(parent_seq)
                bases = np.array(['A', 'C', 'G', 'T'])
                new_seq = np.empty(seq_len, dtype='U1')
                for i, base in enumerate(parent_arr):
                    r = rng.random()
                    if r < p_same:
                        new_seq[i] = base
                    else:
                        alt = bases[bases != base]
                        new_seq[i] = rng.choice(alt)
            else:
                new_seq = parent_seq.copy()
            seq = new_seq.tolist()

        node_to_seq[node] = seq
        for child in node.child_node_iter():
            simulate_node(child, seq)
    simulate_node(tree.seed_node)
    with open(out_fasta, 'w') as f:
        for leaf in tree.leaf_node_iter():
            seq_str = ''.join(node_to_seq[leaf])
            f.write(f">{leaf.taxon.label}\n{seq_str}\n")

    return out_fasta
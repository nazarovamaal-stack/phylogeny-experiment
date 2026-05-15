import os
import sys
import numpy as np
import dendropy
import VGsim
from config import (
    TRANSMISSION_RATE, RECOVERY_RATE, SAMPLING_RATE,
    SUBSTITUTION_WEIGHTS, POPULATION_SIZE,
    EPIDEMIC_TIME, MAX_ITERATIONS, NEUTRAL_LENGTH
)

def run_vgsim_wrapper(n_taxa, n_sites, m, out_fasta, seed, mutation_rate):
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
                adv_haplotype = 'G' * n_sites
                enhanced_rate = TRANSMISSION_RATE * (1.0 + m)
                simulator.set_transmission_rate(enhanced_rate, haplotype=adv_haplotype)
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
                    print(f"Tree file not found, attempt {attempt + 1} failed")
                    continue
            tree = dendropy.Tree.get_from_path(tree_path, schema="newick", rooting="force-unrooted")
        except Exception as e:
            print(f"Simulation/genealogy error: {e}, attempt {attempt+1} failed")
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
            print(f"  Attempt {attempt + 1}: got {len(leaves)} samples (< {n_taxa}), retrying...")
            sys.stdout.flush()
    print(f"Failed to get {n_taxa} samples after {max_attempts} attempts for seed {seed}")
    sys.stdout.flush()
    return None

def run_alisim(tree_path, out_fasta, seed, mutation_rate, seq_len=None):
    if seq_len is None:
        seq_len = NEUTRAL_LENGTH
    tree = dendropy.Tree.get_from_path(tree_path, schema="newick", rooting="force-unrooted")
    for node in tree.postorder_node_iter():
        if node.edge_length is not None:
            node.edge_length = float(node.edge_length)
    for node in tree.postorder_node_iter():
        if node.edge_length is not None:
            node.edge_length *= mutation_rate
    scaled_tree_path = tree_path.replace(".nwk", "_subst.nwk")
    tree.write(path=scaled_tree_path, schema="newick", suppress_rooting=True, unquoted_underscores=True)
    rng = np.random.default_rng(seed)
    node_to_seq = {}

    def simulate_node(node, parent_seq=None):
        if parent_seq is None:
            seq = rng.choice(['A','C','G','T'], size=seq_len)
        else:
            t = node.edge_length if node.edge_length is not None else 0.0
            if t > 0:
                p_same = 0.25 + 0.75 * np.exp(-4.0/3.0 * t)
                new_seq = parent_seq.copy()
                mask = rng.random(seq_len) >= p_same
                if mask.any():
                    bases = np.array(['A','C','G','T'])
                    curr = parent_seq[mask]
                    alt = np.array([bases[bases != c] for c in curr], dtype=object)
                    new_seq[mask] = [rng.choice(a) for a in alt]
                seq = new_seq
            else:
                seq = parent_seq.copy()
        node_to_seq[node] = seq
        for child in node.child_node_iter():
            simulate_node(child, seq)

    simulate_node(tree.seed_node)
    with open(out_fasta, 'w') as f:
        for leaf in tree.leaf_node_iter():
            name = leaf.taxon.label.replace(' ', '_')
            seq_str = ''.join(node_to_seq[leaf])
            f.write(f">{name}\n{seq_str}\n")
    return out_fasta, scaled_tree_path
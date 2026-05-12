"""Филогенетические утилиты"""

import re
import io
import dendropy
from dendropy.calculate import treecompare
from .system import run_cmd
from config import IQTREE_CMD


def run_iqtree(aln_path, out_prefix, seed):
    """
    Запускает IQ-TREE для построения дерева
    """
    run_cmd(
        f"{IQTREE_CMD} -s {aln_path} -m JC -nt AUTO "
        f"-pre {out_prefix} -seed {seed} -redo"
    )
    return f"{out_prefix}.treefile"


def rf_distance(tree1_path, tree2_path, normalized=True):
    """
    Вычисляет Robinson-Foulds расстояние между двумя деревьями.
    normalized=True (по умолчанию) возвращает значение от 0.0 до 1.0.
    """
    taxon_namespace = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get_from_path(
        tree1_path, schema="newick", rooting="force-unrooted",
        taxon_namespace=taxon_namespace
    )
    t2 = dendropy.Tree.get_from_path(
        tree2_path, schema="newick", rooting="force-unrooted",
        taxon_namespace=taxon_namespace
    )
    t2.migrate_taxon_namespace(t1.taxon_namespace)
    rf_raw = treecompare.symmetric_difference(t1, t2)

    if not normalized:
        return rf_raw
    n_taxa = len(t1.taxon_namespace)
    max_rf = 2 * (n_taxa - 3)
    if max_rf <= 0:
        return 0.0

    return rf_raw / max_rf

def branch_score_distance(tree1_path, tree2_path):
    """
    Вычисляет Branch Score Distance (BSD) между двумя деревьями.
    BSD = sqrt( sum_{splits} (len1 - len2)^2 ),
    где для уникальных сплитов длина отсутствующего = 0.
    """
    import re, io, math

    taxon_namespace = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get_from_path(
        tree1_path, schema="newick",
        rooting="force-unrooted", taxon_namespace=taxon_namespace
    )
    with open(tree2_path) as f:
        nw_str = f.read().strip()
    nw_clean = re.sub(r'\[&[^\]]*\]', '', nw_str)
    t2 = dendropy.Tree.get(
        stream=io.StringIO(nw_clean), schema="newick",
        rooting="force-unrooted", taxon_namespace=taxon_namespace
    )
    for leaf in t2.leaf_node_iter():
        if leaf.taxon.label in {t.label for t in t1.taxon_namespace}:
            leaf.taxon = t1.taxon_namespace.get_taxon(leaf.taxon.label)

    def get_splits(tree):
        splits = {}
        for edge in tree.preorder_edge_iter():
            head = edge.head_node
            if head.is_leaf():
                continue
            try:
                leaves = frozenset(leaf.taxon.label for leaf in head.leaf_iter())
            except AttributeError:
                continue
            length = edge.length if edge.length is not None else 0.0
            if leaves not in splits or length != 0.0:
                splits[leaves] = length
        return splits

    splits1 = get_splits(t1)
    splits2 = get_splits(t2)

    all_splits = set(splits1.keys()) | set(splits2.keys())
    bsd_sq = 0.0
    for split in all_splits:
        len1 = splits1.get(split, 0.0)
        len2 = splits2.get(split, 0.0)
        bsd_sq += (len1 - len2) ** 2

    return math.sqrt(bsd_sq)
import re
import io
import dendropy
from dendropy.calculate import treecompare
from .system import run_cmd
from config import IQTREE_CMD

def run_iqtree(aln_path, out_prefix, seed):
    """Запускает IQ-TREE для построения максимально правдоподобного дерева по выравниванию."""
    run_cmd(
        f"{IQTREE_CMD} -s {aln_path} -m JC -nt AUTO "
        f"-pre {out_prefix} -seed {seed} -redo"
    )
    return f"{out_prefix}.treefile"

def rf_distance(tree1_path, tree2_path, normalized=True):
    """Вычисляет расстояние Робинсона-Фолдса между двумя деревьями."""
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

def weighted_rf_distance(tree1_path, tree2_path):
    """Вычисляет взвешенное расстояние Робинсона-Фолдса с учётом длин ветвей."""
    t1 = dendropy.Tree.get_from_path(tree1_path, schema="newick", rooting="force-unrooted")
    t2 = dendropy.Tree.get_from_path(tree2_path, schema="newick", rooting="force-unrooted")
    t1.migrate_taxon_namespace(t2.taxon_namespace)
    return treecompare.weighted_robinson_foulds_distance(t1, t2)
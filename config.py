"""Конфигурация финального эксперимента"""
import os

M_VALUES = [1, 2, 5, 10, 20, 50, 100, 200]
NUM_SELECTED_SITES = [0, 50, 200, 500, 1000]
NEUTRAL_LENGTH = 2000
N_REPLICATES = 5
N_TAXA = 30
WORK_DIR = "exp_baseline"

IQTREE_CMD = "iqtree"
VGSIM_CMD = os.path.join("VGsim", "VGsim_cmd.py")
"""Системные утилиты"""

import os
import subprocess
import sys


def run_cmd(cmd, check=True):
    """
    Выполняет команду в shell или subprocess

    Args:
        cmd: команда (строка или список)
        check: если True, завершает программу при ошибке

    Returns:
        CompletedProcess результат выполнения
    """
    res = subprocess.run(
        cmd,
        shell=True if isinstance(cmd, str) else False,
        capture_output=True,
        text=True
    )

    if check and res.returncode != 0:
        print(f"ошибка: {' '.join(cmd) if isinstance(cmd, list) else cmd}\n{res.stderr}")
        sys.exit(1)

    return res


def ensure_dir(path):
    """Создает директорию, если она не существует"""
    os.makedirs(path, exist_ok=True)
import os
import subprocess
import sys

def run_cmd(cmd, check=True):
    """Выполняет системную команду, возвращает результат. При check=True и ошибке завершает программу."""
    res = subprocess.run(
        cmd,
        shell=True if isinstance(cmd, str) else False,
        capture_output=True,
        text=True
    )
    if check and res.returncode != 0:
        print(f"Error: {cmd}\n{res.stderr}")
        sys.exit(1)
    return res

def ensure_dir(path):
    """Создаёт директорию, если она не существует."""
    os.makedirs(path, exist_ok=True)
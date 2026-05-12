"""Утилиты для работы с последовательностями"""

import os


def parse_fasta(path):
    """
    Парсит FASTA файл в словарь {имя: последовательность}

    Args:
        path: путь к FASTA файлу

    Returns:
        словарь с последовательностями
    """
    seqs = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}

    with open(path) as f:
        name = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            elif name:
                seqs[name].append(line)

    return {k: "".join(v) for k, v in seqs.items()}


def concatenate_fasta(fasta1, fasta2, out_path):
    """
    Конкатенирует два FASTA файла

    Args:
        fasta1: путь к первому FASTA файлу
        fasta2: путь ко второму FASTA файлу
        out_path: путь для выходного файла
    """
    s1 = parse_fasta(fasta1)
    s2 = parse_fasta(fasta2)

    if not s1:
        raise ValueError(f"Пустое нейтральное выравнивание: {fasta1}")

    def normalize_name(name):
        return name.replace(' ', '_')

    s1 = {normalize_name(k): v for k, v in s1.items()}
    s2 = {normalize_name(k): v for k, v in s2.items()}
    with open(out_path, "w") as f:
        for taxon in s1:
            f.write(f">{taxon}\n{s1[taxon]}{s2.get(taxon, '')}\n")
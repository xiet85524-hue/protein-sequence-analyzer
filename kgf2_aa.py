# -*- coding: utf-8 -*-
"""
W1 Day1 · KGF-2 氨基酸组成分析
目标：用纯 Python（不用任何第三方库）读取一段真实的蛋白 FASTA，
     统计每种氨基酸出现了多少次、占百分之几。
运行：python 2026-08-19_W1D1_氨基酸组成.py
"""

# ---------- 1. 读取 FASTA 文件 ----------
def read_fasta(path):
    """逐行读文件，跳过 '>' 开头的标题行，其余行拼成序列字符串。"""
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts)

# 真实 FASTA 路径（用绝对路径，保证从任何位置运行都能找到 Obsidian 里的真实文件）
fasta_path = r"C:\Users\admin\Desktop\实验知识库\08-AI与生物信息\我的蛋白案例库\真实序列\01-KGF-2__一种角质细胞生长因子-透皮肽的融合蛋白在大肠杆菌中的制备与应用_1.fasta"

try:
    seq = read_fasta(fasta_path)
    print("✓ 已从真实 FASTA 文件读取 KGF-2 序列")
except FileNotFoundError:
    # 内置兜底序列（与 FASTA #1 一致），保证换台电脑也能直接跑
    seq = (
        "ALGQDMVSPEATNSSSSSFSSPSSAGRHVRSYNHLQGDVRWRKLFSFTKYFLKIEKNGKVSGTKKENCPY"
        "SILEITSVEIGVVAVKAINSNYYLAMNKKGKLYGSKEFNNDCKLKERIEENGYNTYASFNWQHNGRQMYV"
        "ALNGKGAPRRGQKTRRKNTSAHFLPMVVHSSEQIDN"
    )
    print("⚠ 未找到 FASTA 文件，使用内置 KGF-2 序列")

# ---------- 2. 基础处理 ----------
seq = seq.replace("\n", "").replace(" ", "").upper()
length = len(seq)
print(f"序列长度: {length} 个氨基酸\n")

# ---------- 3. 统计每种氨基酸 ----------
counts = {}
for aa in seq:
    counts[aa] = counts.get(aa, 0) + 1

# ---------- 4. 输出（按数量从多到少） ----------
print(f"{'氨基酸':<6}{'数量':>5}{'百分比':>9}")
print("-" * 22)
for aa, n in sorted(counts.items(), key=lambda x: -x[1]):
    print(f"{aa:<6}{n:>5}{n / length * 100:>8.1f}%")

# ---------- 5. 一致性校验 ----------
total = sum(counts.values())
print("-" * 22)
print(f"计数总和 = {total}，序列长度 = {length}",
      "→ 一致 ✓" if total == length else "→ 不一致 ✗")

# ---------- 6. 给实验员的提示 ----------
if "C" in counts:
    print(f"\n提示：序列含 {counts['C']} 个半胱氨酸(C)，可能形成二硫键，"
          "表达/复性时要注意氧化环境。")
charged = (counts.get("K", 0) + counts.get("R", 0)
           + counts.get("D", 0) + counts.get("E", 0))
print(f"提示：带电荷残基(K/R/D/E)共 {charged} 个，占比 {charged/length*100:.1f}%，"
      "通常利于水溶性。")

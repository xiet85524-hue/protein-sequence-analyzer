# -*- coding: utf-8 -*-
"""
W1 Day2 · 用 Biopython 算 KGF-2 的分子量 / 等电点(pI) / 平均疏水性(GRAVY)
（本文件是 workbuddy-ai-protein 工作副本，FASTA 用绝对路径，双击 run_w1d2.bat 即可运行）
运行：双击同目录 run_w1d2.bat
"""

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

# ---------- 1. 读取 FASTA（绝对路径，保证从任何位置都能找到真实文件） ----------
def read_fasta(path):
    """逐行读文件，跳过 '>' 标题行，其余拼成序列字符串。"""
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts)

fasta_path = r"C:\Users\admin\Desktop\实验知识库\08-AI与生物信息\我的蛋白案例库\真实序列\01-KGF-2__一种角质细胞生长因子-透皮肽的融合蛋白在大肠杆菌中的制备与应用_1.fasta"

try:
    seq = read_fasta(fasta_path)
    print("✓ 已从真实 FASTA 文件读取 KGF-2 序列")
except FileNotFoundError:
    seq = (
        "ALGQDMVSPEATNSSSSSFSSPSSAGRHVRSYNHLQGDVRWRKLFSFTKYFLKIEKNGKVSGTKKENCPY"
        "SILEITSVEIGVVAVKAINSNYYLAMNKKGKLYGSKEFNNDCKLKERIEENGYNTYASFNWQHNGRQMYV"
        "ALNGKGAPRRGQKTRRKNTSAHFLPMVVHSSEQIDN"
    )
    print("⚠ 未找到 FASTA 文件，使用内置 KGF-2 序列")

seq = seq.replace("\n", "").replace(" ", "").upper()

if not BIOPYTHON_OK:
    print("\n✗ 未检测到 Biopython，请先安装： pip install biopython")
    print("（下面只输出能从序列直接算的长度与电荷残基，装好后会显示 MW/pI/GRAVY）\n")
    print(f"序列长度: {len(seq)} aa")
    raise SystemExit(0)

# ---------- 2. 用 Biopython 计算 ----------
a = ProteinAnalysis(seq)

mw = a.molecular_weight()
gravy = a.gravy()
aromaticity = a.aromaticity()
instability = a.instability_index()
net_charge7 = a.charge_at_pH(7.0)

# ---------- 3. 等电点 pI（手写兜底，与 Biopython ProtParam 算法一致） ----------
PKA = {
    "Nterm": 7.5, "Cterm": 3.55,
    "K": 10.0, "R": 12.0, "H": 5.98,
    "D": 4.05, "E": 4.45, "C": 9.0, "Y": 9.84,
}

def net_charge(s, pH):
    c = 0.0
    c += 1.0 / (1.0 + 10 ** (pH - PKA["Nterm"]))
    c += s.count("K") * 1.0 / (1.0 + 10 ** (pH - PKA["K"]))
    c += s.count("R") * 1.0 / (1.0 + 10 ** (pH - PKA["R"]))
    c += s.count("H") * 1.0 / (1.0 + 10 ** (pH - PKA["H"]))
    c -= 1.0 / (1.0 + 10 ** (PKA["Cterm"] - pH))
    c -= s.count("D") * 1.0 / (1.0 + 10 ** (PKA["D"] - pH))
    c -= s.count("E") * 1.0 / (1.0 + 10 ** (PKA["E"] - pH))
    c -= s.count("C") * 1.0 / (1.0 + 10 ** (PKA["C"] - pH))
    c -= s.count("Y") * 1.0 / (1.0 + 10 ** (PKA["Y"] - pH))
    return c

try:
    pI = a.isoelectric_point()
except Exception:
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    pI = (lo + hi) / 2.0

# ---------- 4. 输出 ----------
print("=" * 46)
print("KGF-2 序列特征（Biopython 计算）")
print("=" * 46)
print(f"{'指标':<14}{'数值':>14}{'单位/含义':>18}")
print("-" * 46)
print(f"{'长度':<14}{len(seq):>14}{'aa':>18}")
print(f"{'分子量 MW':<12}{mw/1000:>14.2f}{'kDa':>18}")
print(f"{'等电点 pI':<12}{pI:>14.2f}{'(净电荷=0的pH)':>18}")
print(f"{'平均疏水性GRAVY':<10}{gravy:>14.3f}{'负=亲水':>18}")
print(f"{'芳香族比例':<12}{aromaticity:>14.3f}{'(F/Y/W)':>18}")
print(f"{'不稳定指数':<12}{instability:>14.1f}{'(<40较稳定)':>18}")
print(f"{'净电荷@pH7':<12}{net_charge7:>14.2f}{'(带电与否)':>18}")
print("-" * 46)

# ---------- 5. 给实验员的提示 ----------
charged = seq.count("K") + seq.count("R") + seq.count("D") + seq.count("E")
cys = seq.count("C")
print("\n→ 实验端怎么用这三个数：")
print(f"  · 分子量 {mw/1000:.1f} kDa：SDS-PAGE 上大概在这个位置出条带；"
      "电泳/超滤膜截留选型可参考。")
print(f"  · pI = {pI:.2f}：用 Q/SP 离子交换时，缓冲液 pH 应离 pI 至少 1 个单位"
      f"（例如 pH {pI+1.5:.1f} 或 {pI-1.5:.1f}）让蛋白带净电荷、好结合；"
      "也提示最易聚集沉淀的 pH 就在 pI 附近。")
if gravy < 0:
    print(f"  · GRAVY = {gravy:.3f}（负值=亲水）：整体偏亲水，对可溶性是好事；"
          "但仍有疏水区需注意。")
else:
    print(f"  · GRAVY = {gravy:.3f}（正值=偏疏水）：有一定包涵体/聚集倾向，"
          "表达时留意可溶性。")
print(f"  · 带电荷残基(K/R/D/E)共 {charged} 个 → 水溶性通常较好。")
if cys >= 2:
    print(f"  · 含 {cys} 个 Cys → 可能形成二硫键或错误配对导致聚集，"
          "复性/纯化注意氧化环境控制。")
if instability < 40:
    print(f"  · 不稳定指数 {instability:.1f} < 40 → 通常较稳定。")
else:
    print(f"  · 不稳定指数 {instability:.1f} ≥ 40 → 提示一定不稳定性风险，"
          "可结合后面 W3–W4 的结构预测一起看。")

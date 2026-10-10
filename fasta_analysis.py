# -*- coding: utf-8 -*-
"""
W1 Day4 · 可复用蛋白序列分析工具（升级版）
===========================================
在 W1 Day3 基础上：
  1) 把 KNOWN 扩到「案例库里全部真实 FASTA」——一次循环批量出全部对比表。
  2) 新增 batch()：循环分析所有已知蛋白，打印横排对比，并导出 CSV（科研资产）。
  3) 新增「序列类型自动判别」：含大量 A/T/G/C（疑似核酸）或含非标准字母的序列
     会标记为「疑似核酸/非蛋白」，跳过错误的蛋白指标统计，避免误判（KGK-2 变体3 就是教训）。
  4) 新增「可扩展指标区」：示范加一个对你表达/复性最有用新指标
     （潜在二硫键对数 + 疏水核心残基占比），其余留给你自己补（见 extra_metrics()）。

用法（双击对应 bat 即可，也可命令行）：
  python fasta_analysis.py            # 默认分析 KGF-2（变体1）
  python fasta_analysis.py opn        # 分析骨桥蛋白 OPN（变体1）
  python fasta_analysis.py kegf       # 分析 KEGF2-5+V123I（变体1）
  python fasta_analysis.py all        # ← W1 Day4 新增：批量分析全部真实序列，并导出 CSV
  python fasta_analysis.py "某路径.fasta"   # 高级：直接给任意 FASTA 路径

输出：氨基酸组成 + 分子量/等电点/疏水性/不稳定指数/净电荷 + 实验端提示
"""

import sys
import os
import csv

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

# ---------- 已知蛋白的 FASTA 路径（覆盖案例库全部真实序列，第三元组=类型） ----------
BASE = r"C:\Users\admin\Desktop\实验知识库\08-AI与生物信息\我的蛋白案例库\真实序列"
# KNOWN[code] = (path, name, kind)   kind: "prot"=蛋白  "na"=疑似核酸（跳过蛋白统计）
KNOWN = {
    # —— KGF-2 三个变体（原 kgf2 = 变体1）——
    "kgf2_1": (os.path.join(BASE, "01-KGF-2__一种角质细胞生长因子-透皮肽的融合蛋白在大肠杆菌中的制备与应用_1.fasta"), "KGF-2（变体1·176aa）", "prot"),
    "kgf2_2": (os.path.join(BASE, "01-KGF-2__一种角质细胞生长因子-透皮肽的融合蛋白在大肠杆菌中的制备与应用_2.fasta"), "KGF-2（变体2·197aa）", "prot"),
    "kgf2_3": (os.path.join(BASE, "01-KGF-2__一种角质细胞生长因子-透皮肽的融合蛋白在大肠杆菌中的制备与应用_3.fasta"), "KGF-2（变体3·DNA序列!）", "na"),
    # —— 骨桥蛋白 OPN 三个变体（原 opn = 变体1）——
    "opn_1":  (os.path.join(BASE, "03-骨桥蛋白OPN__AI辅助设计-骨桥蛋白序列_1.fasta"), "骨桥蛋白 OPN（变体1·169aa）", "prot"),
    "opn_2":  (os.path.join(BASE, "03-骨桥蛋白OPN__AI辅助设计-骨桥蛋白序列_2.fasta"), "骨桥蛋白 OPN（变体2·314aa）", "prot"),
    "opn_3":  (os.path.join(BASE, "03-骨桥蛋白OPN__AI辅助设计-骨桥蛋白序列_3.fasta"), "骨桥蛋白 OPN（变体3·278aa）", "prot"),
    # —— KEGF2-5+V123I 两个变体（原 kegf = 变体1）——
    "kegf_1": (os.path.join(BASE, "05-生长因子系列__KEGF2-5+V123I序列_1.fasta"), "KEGF2-5+V123I（变体1·172aa）", "prot"),
    "kegf_2": (os.path.join(BASE, "05-生长因子系列__KEGF2-5+V123I序列_2.fasta"), "KEGF2-5+V123I（变体2·173aa）", "prot"),
    # —— 分子生物学：Cathepsin B 酶切序列 ——
    "catb":   (os.path.join(BASE, "07-分子生物学__Cathepsin B酶切序列_1.fasta"), "Cathepsin B 酶切序列（171aa）", "prot"),
}

# 兼容旧代号（Day3 的 bat 仍在用）
ALIAS = {"kgf2": "kgf2_1", "opn": "opn_1", "kegf": "kegf_1"}

# KGF-2 内置兜底序列（找不到文件时用）
FALLBACK = (
    "ALGQDMVSPEATNSSSSSFSSPSSAGRHVRSYNHLQGDVRWRKLFSFTKYFLKIEKNGKVSGTKKENCPY"
    "SILEITSVEIGVVAVKAINSNYYLAMNKKGKLYGSKEFNNDCKLKERIEENGYNTYASFNWQHNGRQMYV"
    "ALNGKGAPRRGQKTRRKNTSAHFLPMVVHSSEQIDN"
)

# 疏水核心残基（经典疏水簇，易形成疏水斑/聚集核）
HYDROPHOBIC = set("AVILMFWYC")

# 非标准氨基酸字母（蛋白里不该大量出现；X=未知，U=硒代半胱氨酸极少见）
NON_AA = set("BJOUXZ")


# ---------- 序列类型自动判别（★ W1 Day4 新增：数据清洗防线） ----------
def looks_like_protein(seq):
    """粗略判断序列是否像蛋白：若 A/T/G/C 占比过高→疑似核酸；含非标准字母→非蛋白。"""
    s = seq.upper()
    if not s:
        return False
    nt = sum(s.count(c) for c in "ATGC")
    if nt / len(s) > 0.9:          # 真蛋白里只有 A(丙氨酸)属于 ATGC，占比通常很低
        return False
    if any(ch in NON_AA for ch in s):
        return False
    return True


# ---------- 1. 读取 FASTA ----------
def read_fasta(path):
    """逐行读文件，跳过 '>' 标题行，其余拼成序列字符串。"""
    parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts.append(line)
    return "".join(parts)


# ---------- 2. pI 兜底计算（手写，与 Biopython 算法一致） ----------
PKA = {
    "Nterm": 7.5, "Cterm": 3.55,
    "K": 10.0, "R": 12.0, "H": 5.98,
    "D": 4.05, "E": 4.45, "C": 9.0, "Y": 9.84,
}

def _net_charge(s, pH):
    c = 1.0 / (1.0 + 10 ** (pH - PKA["Nterm"]))
    c += s.count("K") * 1.0 / (1.0 + 10 ** (pH - PKA["K"]))
    c += s.count("R") * 1.0 / (1.0 + 10 ** (pH - PKA["R"]))
    c += s.count("H") * 1.0 / (1.0 + 10 ** (pH - PKA["H"]))
    c -= 1.0 / (1.0 + 10 ** (PKA["Cterm"] - pH))
    c -= s.count("D") * 1.0 / (1.0 + 10 ** (PKA["D"] - pH))
    c -= s.count("E") * 1.0 / (1.0 + 10 ** (PKA["E"] - pH))
    c -= s.count("C") * 1.0 / (1.0 + 10 ** (PKA["C"] - pH))
    c -= s.count("Y") * 1.0 / (1.0 + 10 ** (PKA["Y"] - pH))
    return c

def _pI_fallback(s):
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if _net_charge(s, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


# ========== 可扩展指标区（★ W1 Day4：你亲手加指标的地方） ==========
# analyze() 已经算好了长度/组成/MW/pI/GRAVY/芳香族/不稳定/净电荷/带电荷/Cys。
# 下面示范加一个「对你表达/复性最有用」的新指标：
#   · cys_pairs：潜在二硫键对数 = Cys 数 // 2（≥2 提示氧化折叠/复性要控环境）
#   · hydro_ratio：疏水核心残基占比（越高越易形成疏水斑/包涵体倾向）
# 然后 extra_metrics() 是一个「留给你自己实现」的钩子——你可以加任意新指标，
# 比如：带正电残基占比、脯氨酸（柔韧/转角）占比、N/Q 脱酰胺风险位点计数……
# 只要 return 一个 dict，batch() 会自动把它拼进对比表最后一列。
def extra_metrics(seq):
    """★ 占位：你亲手加的新指标。返回 dict（键=指标名，值=数字或字符串）。"""
    s = seq.upper()
    n_cys = s.count("C")
    hydro = sum(s.count(a) for a in HYDROPHOBIC)
    out = {
        "eps280": int(round(s.count("W") * 5500 + s.count("Y") * 1490 + (s.count("C") // 2) * 125)),
        "潜在二硫键对数": n_cys // 2,
        "疏水核心残基占比(%)": round(hydro / len(s) * 100, 1),
        # ↓↓↓ 取消注释，或自己加一行，就是你亲手加的指标 ↓↓↓
        # "脯氨酸占比(%)": round(s.count("P") / len(s) * 100, 1),
    }
    return out
# ========================= 可扩展指标区结束 =========================


# ---------- 3. 核心：分析一条序列 ----------
def analyze(seq):
    """返回所有指标（dict）。这是你第一个"自己的分析函数"。"""
    seq = seq.replace("\n", "").replace(" ", "").upper()
    is_prot = looks_like_protein(seq)
    res = {"length": len(seq), "is_protein": is_prot}
    if not is_prot:
        res["warning"] = "疑似核酸/非蛋白序列，已跳过蛋白指标统计（请核对原始数据）"
        return res

    # 氨基酸组成（按数量降序）
    comp = {}
    for aa in seq:
        comp[aa] = comp.get(aa, 0) + 1
    res["composition"] = dict(sorted(comp.items(), key=lambda x: -x[1]))

    if not BIOPYTHON_OK:
        res["biopython"] = False
        res.update(extra_metrics(seq))
        return res

    res["biopython"] = True
    a = ProteinAnalysis(seq)
    res["mw"] = a.molecular_weight()
    res["gravy"] = a.gravy()
    res["aromaticity"] = a.aromaticity()
    res["instability"] = a.instability_index()
    res["net_charge7"] = a.charge_at_pH(7.0)
    try:
        res["pI"] = a.isoelectric_point()
    except Exception:
        res["pI"] = _pI_fallback(seq)
    res["charged"] = seq.count("K") + seq.count("R") + seq.count("D") + seq.count("E")
    res["cys"] = seq.count("C")
    # 叠加你亲手加（或示范）的新指标
    res.update(extra_metrics(seq))
    return res


# ---------- 4. 打印单个报告 ----------
def report(seq, name="未知蛋白"):
    r = analyze(seq)
    print("=" * 50)
    print(f"蛋白：{name}")
    print("=" * 50)
    if not r.get("is_protein", True):
        print(f"[WARN] {r.get('warning','')}")
        print(f"  序列长度 = {r['length']}（按字符计，可能为核苷酸数）")
        print("  → 建议：回到原始 docx/数据库核对，这条序列不该进蛋白分析。")
        return r
    print(f"{'指标':<14}{'数值':>14}{'单位/含义':>20}")
    print("-" * 50)
    print(f"{'长度':<14}{r['length']:>14}{'aa':>20}")
    if r.get("biopython"):
        print(f"{'分子量 MW':<12}{r['mw']/1000:>14.2f}{'kDa':>20}")
        print(f"{'等电点 pI':<12}{r['pI']:>14.2f}{'(净电荷=0的pH)':>20}")
        print(f"{'平均疏水性GRAVY':<10}{r['gravy']:>14.3f}{'负=亲水':>20}")
        print(f"{'芳香族比例':<12}{r['aromaticity']:>14.3f}{'(F/Y/W)':>20}")
        print(f"{'不稳定指数':<12}{r['instability']:>14.1f}{'(<40较稳定)':>20}")
        print(f"{'净电荷@pH7':<12}{r['net_charge7']:>14.2f}{'(带电与否)':>20}")
        print(f"{'带电荷残基':<12}{r['charged']:>14}{'(K/R/D/E)':>20}")
        print(f"{'Cysteine':<14}{r['cys']:>14}{'二硫键风险':>20}")
    else:
        print("[WARN] 未检测到 Biopython，仅输出长度与组成；装好 biopython 后会显示 MW/pI/GRAVY。")
    print("-" * 50)

    print("\n→ 氨基酸组成（按数量降序）：")
    for aa, n in r["composition"].items():
        print(f"  {aa}: {n:>3}  ({n/r['length']*100:5.1f}%)")

    if r.get("biopython"):
        print("\n→ 实验端怎么用这些数：")
        print(f"  · 分子量 {r['mw']/1000:.1f} kDa：SDS-PAGE 条带位置参考；超滤膜/层析柱选型参考。")
        print(f"  · pI = {r['pI']:.2f}：Q/SP 离子交换时缓冲液 pH 应离 pI 至少 1 个单位"
              f"（例如 pH {r['pI']+1.5:.1f} 或 {r['pI']-1.5:.1f}），也提示最易聚集的 pH。")
        if r["gravy"] < 0:
            print(f"  · GRAVY = {r['gravy']:.3f}（亲水）：对可溶性是好事，但仍有疏水区需注意。")
        else:
            print(f"  · GRAVY = {r['gravy']:.3f}（偏疏水）：有一定包涵体/聚集倾向，表达时留意可溶性。")
        print(f"  · 带电荷残基(K/R/D/E)共 {r['charged']} 个 → 水溶性通常较好。")
        if r["cys"] >= 2:
            print(f"  · 含 {r['cys']} 个 Cys → 可能形成二硫键或错误配对导致聚集，复性/纯化注意氧化环境。")
        if r["instability"] < 40:
            print(f"  · 不稳定指数 {r['instability']:.1f} < 40 → 通常较稳定。")
        else:
            print(f"  · 不稳定指数 {r['instability']:.1f} ≥ 40 → 提示一定不稳定性风险，可结合结构预测一起看。")
    return r


# ---------- 5. 批量：循环分析全部已知蛋白 + 导出 CSV ----------
def batch():
    print("=" * 70)
    print("[INFO] 批量分析：案例库全部真实序列特征对比")
    print("=" * 70)
    rows = []
    for code, (path, name, kind) in KNOWN.items():
        try:
            seq = read_fasta(path)
        except FileNotFoundError:
            seq = FALLBACK
            name = name + "（未找到文件，用内置序列）"
        r = analyze(seq)
        if not r.get("is_protein", True):
            print(f"\n[{code}] {name}")
            print(f"  [WARN] {r.get('warning','')}（长度={r['length']}）")
            print("     → 已跳过，不计入蛋白对比表（请回原始数据核对）。")
            continue
        rows.append((code, name, r))
        print(f"\n[{code}] {name}")
        print(f"  长度={r['length']}aa | MW={r.get('mw',0)/1000:.2f}kDa | "
              f"pI={r.get('pI',0):.2f} | GRAVY={r.get('gravy',0):.3f} | "
              f"不稳定={r.get('instability',0):.1f} | Cys={r.get('cys',0)} | "
              f"二硫键对={r.get('潜在二硫键对数','-')} | 疏水%={r.get('疏水核心残基占比(%)','-')}")

    # 导出 CSV（科研资产：一张可沉淀的对比数据表）
    print("[SCRIPT_DIR]", os.path.dirname(os.path.abspath(__file__)))
    out_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "w1d4_序列特征对比表.csv")
    fields = ["代号", "名称", "长度", "分子量kDa", "pI", "GRAVY", "芳香族",
              "不稳定指数", "净电荷@pH7", "带电荷残基", "Cys",
              "潜在二硫键对数", "疏水核心残基占比(%)", "eps280"]
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for code, name, r in rows:
            w.writerow([
                code, name, r["length"],
                round(r.get("mw", 0) / 1000, 2), round(r.get("pI", 0), 2),
                round(r.get("gravy", 0), 3), round(r.get("aromaticity", 0), 3),
                round(r.get("instability", 0), 1), round(r.get("net_charge7", 0), 2),
                r.get("charged", 0), r.get("cys", 0),
                r.get("潜在二硫键对数", "-"), r.get("疏水核心残基占比(%)", "-"), r.get("eps280", "-"),
            ])
    print("\n" + "=" * 70)
    print(f"[OK] 蛋白对比表已导出：{out_csv}（共 {len(rows)} 条蛋白，已自动剔除疑似核酸序列）")
    print("（这就是今天沉淀的科研资产：一份可直接进 Excel/案例库的数据表）")
    return rows


# ---------- 6. 主程序 ----------
def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else "kgf2"
    # 别名兼容（kgf2/opn/kegf）
    key = ALIAS.get(arg, arg)

    if key == "all":
        batch()
        return

    if key in KNOWN:
        path, name, kind = KNOWN[key]
        try:
            seq = read_fasta(path)
            print(f"[OK] 已从真实 FASTA 读取：{name}")
        except FileNotFoundError:
            seq = FALLBACK
            name = name + "（未找到文件，用内置序列）"
            print("[WARN] 未找到 FASTA 文件，使用内置 KGF-2 序列")
    else:
        try:
            seq = read_fasta(arg)
            name = os.path.basename(arg)
            print(f"[OK] 已从文件读取：{name}")
        except FileNotFoundError:
            print(f"[ERROR] 找不到文件或未知蛋白代号：{arg}")
            print("可用代号：", ", ".join(list(KNOWN.keys()) + ["all"]))
            return
    report(seq, name=name)


if __name__ == "__main__":
    main()




def count_nxs(seq):
    n = 0
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
            n = n + 1
    return n
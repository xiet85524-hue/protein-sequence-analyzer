import os
import fasta_analysis as fa

MASS = {"A":71.0788,"R":156.1875,"N":114.1038,"D":115.0886,"C":103.1388,
        "E":129.1155,"Q":128.1307,"G":57.0519,"H":137.1411,"I":113.1594,
        "L":113.1594,"K":128.1741,"M":131.1926,"F":147.1766,"P":97.1167,
        "S":87.0782,"T":101.1051,"W":186.2132,"Y":163.1760,"V":99.1326}

def mw(seq):
    total = 18.01528
    for a in seq:
        total = total + MASS.get(a, 0.0)
    return total

mw_list = []
for name in sorted(os.listdir(fa.BASE)):
    if not name.endswith(".fasta"):
        continue
    seq = fa.read_fasta(os.path.join(fa.BASE, name)).upper()
    if not fa.looks_like_protein(seq):
        continue
    mw_list.append(round(mw(seq) / 1000, 2))

print("[COUNT]", len(mw_list))
mw_list.sort(reverse=True)
for m in mw_list:
    print("[MW-kDa]", m)

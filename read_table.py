import pandas as pd

df = pd.read_csv("w1d4_序列特征对比表.csv", encoding="utf-8-sig")
mask = df["长度"] > 200
print("开关为 True 的个数 =", mask.sum())
big = df[mask]
print(big[["代号", "长度"]])

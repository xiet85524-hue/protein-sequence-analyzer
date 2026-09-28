# 蛋白序列分析工具（Protein Sequence Analyzer）

一句话：把一条蛋白 FASTA 序列丢进去，输出它的 ___ 和 ___

## 解决什么问题
以前我要算一个蛋白的分子量 / 等电点，得 ___（手算？查网页？一个个点？）……
这个脚本把它变成了 ___

## 怎么运行
需要 Python 3.13 + biopython。

    python fasta_analysis.py          # 默认分析 KGF-2（变体1，176 aa）
    python fasta_analysis.py opn      # 分析骨桥蛋白 OPN
    python fasta_analysis.py all      # 批量分析案例库全部真实序列，并导出对比表
    python fasta_analysis.py "某路径.fasta"   # 任意 FASTA 文件

## 输出什么
屏幕上：氨基酸组成 + 分子量 / 等电点 / 疏水性(GRAVY) / 不稳定指数 / 净电荷 / 潜在二硫键 /
        N-X-S/T 糖基化位点数 + 实验端提示
文件里：w1d4_序列特征对比表.csv（UTF-8-SIG，可直接用 Excel 打开）

## 用了什么数据
我的真实蛋白案例库，共 ___ 条 FASTA。
其中 ___ 那条其实是一段 DNA —— 脚本会自动识别并跳过（这是我踩过的真实坑）。

## 我自己真正理解的部分 / 还没搞懂的部分
- 懂：___
- 还没搞懂：___

## 已知局限（诚实写出来）
- `run_w1d4_all.bat` 里的路径是写死的绝对路径，换台电脑就用不了
- ___

## 待办
- [ ] 补第 9、10 个指标：___
- [ ] 推到 GitHub
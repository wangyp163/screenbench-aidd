#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 0 · random split vs scaffold split 泄漏实测

问题：同样的模型同样的特征，随机切分 vs 骨架切分（scaffold split）的 ROC-AUC 差多少？
差值 = 「测试集和训练集长得像」带来的泄漏虚高——真实场景遇到的是新骨架，scaffold 分数才可信。

数据：MoleculeNet HIV（TDC 基准同源，deepchemdata S3 下载，SMILES 去重后 ~41k）
特征：Morgan ECFP4 (r=2, 2048 bit)；模型：RandomForest(200 trees)
运行：python phase0_split_leakage.py  （CPU，约 10 分钟）
"""
import os
import time
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

RDLogger.DisableLog("rdApp.*")
N_SEEDS = 3
N_TREES = 200
N_JOBS = 14

print("加载数据（MoleculeNet HIV · TDC 同源）...")
df = pd.read_csv("data/HIV.csv")
n_raw = len(df)
df = df.drop_duplicates(subset=["smiles"]).reset_index(drop=True)
print(f"原始 {n_raw} 行，按 SMILES 去重后 {len(df)} 行，活性占比 {df['HIV_active'].mean():.1%}")

fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps, scaffolds, labels = [], [], []
bad = 0
for smi, yv in zip(df["smiles"], df["HIV_active"]):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        bad += 1
        continue
    fps.append(np.array(fp_gen.GetFingerprint(mol), dtype=np.uint8))
    scaffolds.append(MurckoScaffold.MurckoScaffoldSmiles(mol=mol))
    labels.append(int(yv))
X = np.asarray(fps, dtype=np.float32)
y = np.asarray(labels)
groups = pd.Series(scaffolds)
n_scaffold = groups.nunique()
print(f"解析成功 {len(y)} 条（失败 {bad}），活性占比 {y.mean():.1%}，唯一骨架 {n_scaffold} 个")

def random_auc(seed):
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, stratify=y, random_state=seed)
    clf = RandomForestClassifier(n_estimators=N_TREES, n_jobs=N_JOBS, random_state=seed)
    clf.fit(Xtr, ytr)
    return roc_auc_score(yte, clf.predict_proba(Xte)[:, 1])

def scaffold_auc(seed):
    # 按骨架分组切分：同一骨架的分子只出现在训练或测试一侧
    rng = np.random.RandomState(seed)
    uniq = pd.unique(groups)
    rng.shuffle(uniq)
    hold = set(uniq[: int(len(uniq) * 0.2)])
    te = np.isin(groups.values, list(hold))
    clf = RandomForestClassifier(n_estimators=N_TREES, n_jobs=N_JOBS, random_state=seed)
    clf.fit(X[~te], y[~te])
    assert y[te].sum() > 0, "测试集无正样本"
    return roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1])

results = {"random": [], "scaffold": []}
for seed in range(N_SEEDS):
    for name, fn in (("random", random_auc), ("scaffold", scaffold_auc)):
        t0 = time.time()
        auc = fn(seed)
        results[name].append(auc)
        print(f"seed={seed} {name:9s} ROC-AUC = {auc:.4f}   ({time.time()-t0:.0f}s)")

r = np.array(results["random"])
s = np.array(results["scaffold"])
gap = r.mean() - s.mean()
report = f"""# Phase 0 · 切分方式泄漏实测（MoleculeNet HIV · TDC 同源）

- 数据：{len(y)} 分子（活性 {y.mean():.1%}），唯一骨架 {n_scaffold} 个
- 特征：Morgan ECFP4 (2048 bit)；模型：RandomForest({N_TREES} trees)；{N_SEEDS} 个种子
- 测试集均 ~20%：random = 随机分层；scaffold = 按骨架分组（同骨架不跨侧）

| 切分 | ROC-AUC (mean ± std) |
|---|---|
| random   | {r.mean():.4f} ± {r.std():.4f} |
| scaffold | {s.mean():.4f} ± {s.std():.4f} |
| **泄漏虚高** | **+{gap:.4f}（相对 {gap/s.mean():.1%}）** |

## 一句话结论
随机切分比骨架切分虚高 **{gap/s.mean():.1%}**——高出的部分不是模型学到了通用规律，
而是「测试分子和训练分子同骨架」的记忆分。真实筛选场景遇到的是新骨架，所以药物域统一评测
（TDC 等）默认 scaffold split。这就是「冻结评测集 + 防泄漏切分」存在的意义。
"""
print("\n" + report)
os.makedirs("results", exist_ok=True)
with open("results/phase0.md", "w") as f:
    f.write(report)
print("已写入 results/phase0.md")

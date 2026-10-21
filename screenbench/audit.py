"""数据审计：任意分子 CSV 的骨架多样性 + 集内近重复检测（Tanimoto > 0.9）。"""
import sys
from collections import Counter

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


def audit_csv(path: str, smiles_col: str = "smiles") -> str:
    df = pd.read_csv(path)
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps, scaffolds = [], []
    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        fps.append(np.array(fp_gen.GetFingerprint(mol), dtype=np.uint8))
        scaffolds.append(MurckoScaffold.MurckoScaffoldSmiles(mol=mol))
    X = np.asarray(fps, dtype=bool)
    n = len(X)
    top_scaff = Counter(scaffolds).most_common(5)

    # 近重复：全对 O(n^2) 矩阵乘法，>2 万分子时抽样到 2 万
    if n > 20000:
        idx = np.random.RandomState(0).choice(n, 20000, replace=False)
        X = X[idx]
        n = 20000
    dots = (X @ X.T).astype(np.float32)
    tan = dots / (2048 * 2 - dots + 1e-9)
    iu = np.triu_indices(n, k=1)
    near = float((tan[iu] > 0.9).mean())

    lines = [f"# audit · {path}", "",
             f"- 有效分子 {len(scaffolds)}，唯一骨架 {len(set(scaffolds))}（{len(set(scaffolds))/max(1,len(scaffolds)):.1%}）",
             "- Top5 骨架：" + ", ".join(f"{s}({c})" for s, c in top_scaff),
             f"- 集内近重复对（Tanimoto>0.9）占比：{near:.4%}",
             "", "> 近重复占比高 = 数据集信息量虚高，任何 split 的分数都要打折看待。"]
    return "\n".join(lines)

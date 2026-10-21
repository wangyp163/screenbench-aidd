"""Redocking 验证：晶体配体重新对接回自己的口袋，RMSD < 2 Å 判定管线正确性。

原理：盒子就是用这个配体定义的（self-docking），如果制备链路（obabel 受体 + Meeko 配体 +
盒子规则 + vina）有系统性错误，best pose 不可能回到晶体位置。这是 VS 管线的"单元测试"。

RMSD 计算不绕道 obabel/sdf：vina 输出 pdbqt 保持输入配体的原子顺序，
直接对 input/out 的 MODEL 1 重原子坐标逐对求 RMSD。
"""
import json
import math
import time
from pathlib import Path

from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from . import prep, engines, registry
from .paths import REPO

RUNS = REPO / "runs" / "redock"


def _heavy_coords(path: Path, model: int = 1) -> list[tuple[float, float, float]]:
    coords, in_model, saw_model = [], False, False
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            in_model = int(line.split()[1]) == model
            continue
        if line.startswith("ENDMDL"):
            in_model = False
            continue
        if saw_model and not in_model:
            continue
        if line.startswith(("ATOM", "HETATM")):
            atype = line[77:79].strip()
            if atype and not atype.startswith("H"):
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


def _rmsd(a, b) -> float:
    assert len(a) == len(b), f"原子数不匹配 {len(a)} vs {len(b)}（拓扑被改动？）"
    return math.sqrt(sum((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                         for (x1, y1, z1), (x2, y2, z2) in zip(a, b)) / len(a))


def redock_and_report(set_name: str, tg: str, exhaustiveness: int = 8, threshold: float = 2.0):
    tg_dir = registry.DEKOIS_ROOT / tg
    cache = REPO / "cache" / set_name / tg / "redock"
    cache.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # 1. 晶体配体 → pdbqt（坐标即晶体 pose）
    lig_pdbqt = cache / "xtal_lig.pdbqt"
    if not lig_pdbqt.is_file():
        mol = Chem.SDMolSupplier(str(tg_dir / "ligand.sdf"), sanitize=True, removeHs=False)[0]
        if any(a.GetTotalNumHs() for a in mol.GetAtoms()):
            mol = Chem.AddHs(mol, addCoords=True)
        pdbqt, ok, _ = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
        assert ok, "meeko 转换失败"
        lig_pdbqt.write_text(pdbqt)
    ref = _heavy_coords(lig_pdbqt)

    # 2. 盒子 + docking（单配体，高 exhaustiveness）
    box = prep.box_from_ligand_sdf(tg_dir)
    rec = prep.prep_receptor_dekois(tg_dir, cache.parent)
    out_dir = cache / f"dock_exh{exhaustiveness}"
    lig_dir = cache / "lig"
    lig_dir.mkdir(exist_ok=True)
    (lig_dir / "xtal_lig.pdbqt").write_text(lig_pdbqt.read_text())
    engines.run_vina(rec, lig_dir, out_dir, box, seed=42, exhaustiveness=exhaustiveness, n_parallel=1)

    # 3. best pose RMSD
    out_file = out_dir / "xtal_lig_out.pdbqt"
    probe = _heavy_coords(out_file, model=1)
    scores = engines.parse_results(out_dir)
    result = {
        "set": set_name, "target": tg, "engine": engines.engine_version(),
        "exhaustiveness": exhaustiveness, "n_heavy_atoms": len(ref),
        "best_affinity": scores.get("xtal_lig"),
        "rmsd_A": round(_rmsd(ref, probe), 3),
        "threshold_A": threshold,
        "pass": _rmsd(ref, probe) < threshold,
        "seconds": round(time.time() - t0, 1),
    }
    RUNS.mkdir(parents=True, exist_ok=True)
    (RUNS / f"{set_name}_{tg}.json").write_text(json.dumps(result, indent=2))
    print(f"[redock] {tg}: RMSD={result['rmsd_A']}Å (<{threshold}Å) → {'PASS ✓' if result['pass'] else 'FAIL ✗'}"
          f"  affinity={result['best_affinity']}  用时 {result['seconds']}s")
    return result

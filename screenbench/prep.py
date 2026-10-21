"""前处理：受体 pdbqt（obabel）与配体 pdbqt（RDKit+Meeko），带缓存。"""
import gzip
import re
import subprocess
from pathlib import Path

from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from .paths import tool_bin, tool_env, REPO
from . import registry


def prep_receptor(set_name: str, tg: str, cache_dir: Path) -> Path:
    out = cache_dir / "rec.pdbqt"
    if out.is_file():
        return out
    cache_dir.mkdir(parents=True, exist_ok=True)
    src = registry.DATA / set_name / tg / f"{tg}_pdb.pdb.gz"
    raw = gzip.open(src, "rt", errors="ignore").read()
    # 只保留蛋白主链记录：去水/去杂环配体/去 CONECT（DUD-E 受体 pdb 偶带 HETATM）
    clean = "\n".join(l for l in raw.splitlines()
                      if l.startswith(("ATOM", "TER", "END", "HEADER", "TITLE")))
    pdb = cache_dir / "rec_clean.pdb"
    pdb.write_text(clean + "\n")
    r = subprocess.run(
        [tool_bin("obabel"), str(pdb), "-O", str(out), "-xr"],
        capture_output=True, text=True, env=tool_env())
    if not out.is_file():
        raise RuntimeError(f"obabel 受体制备失败: {r.stdout}\n{r.stderr}")
    return out


def box_from_crystal_ligand(set_name: str, tg: str, pad_A: float = 4.0, max_edge_A: float = 26.0) -> dict:
    """盒子规则（记录进 run config，保证可复现）：晶体配体 bbox + pad，封顶边长。"""
    src = registry.DATA / set_name / tg / f"{tg}_lig.mol2.gz"
    mol = None
    for _, block in registry.read_mol2_gz(src):
        m = Chem.MolFromMol2Block(block, sanitize=True, removeHs=False)
        if m is not None and m.GetNumConformers() > 0:
            mol = m
            break
    if mol is None:
        raise RuntimeError(f"{tg} 晶体配体解析失败，无法定义盒子")
    conf = mol.GetConformer()
    pts = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    xs, ys, zs = zip(*pts)

    def axis(lo, hi):
        size = min((hi - lo) + 2 * pad_A, max_edge_A)
        center = (hi + lo) / 2
        return round(center, 3), round(size, 3)

    cx, sx = axis(min(xs), max(xs))
    cy, sy = axis(min(ys), max(ys))
    cz, sz = axis(min(zs), max(zs))
    return {"rule": "crystal_ligand_bbox", "center_x": cx, "center_y": cy, "center_z": cz,
            "size_x": sx, "size_y": sy, "size_z": sz, "pad_A": pad_A, "max_edge_A": max_edge_A}


_ID_RE = re.compile(r"^[A-Za-z0-9_\-]+$")


# ── DEKOIS(airdd) 布局：receptor.pdb / ligand.sdf / actives.sdf / decoys.sdfgz ──

def prep_receptor_dekois(tg_dir: Path, cache_dir: Path) -> Path:
    out = cache_dir / "rec.pdbqt"
    if out.is_file():
        return out
    cache_dir.mkdir(parents=True, exist_ok=True)
    r = subprocess.run([tool_bin("obabel"), str(tg_dir / "receptor.pdb"), "-O", str(out), "-xr"],
                       capture_output=True, text=True, env=tool_env())
    if not out.is_file():
        raise RuntimeError(f"obabel 受体制备失败: {r.stdout}\n{r.stderr}")
    return out


def box_from_ligand_sdf(tg_dir: Path, pad_A: float = 4.0, max_edge_A: float = 26.0) -> dict:
    mol = Chem.SDMolSupplier(str(tg_dir / "ligand.sdf"), sanitize=True, removeHs=False)[0]
    if mol is None or mol.GetNumConformers() == 0:
        raise RuntimeError("晶体配体解析失败")
    conf = mol.GetConformer()
    xs, ys, zs = zip(*[list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

    def axis(lo, hi):
        return round((hi + lo) / 2, 3), round(min((hi - lo) + 2 * pad_A, max_edge_A), 3)

    cx, sx = axis(min(xs), max(xs))
    cy, sy = axis(min(ys), max(ys))
    cz, sz = axis(min(zs), max(zs))
    return {"rule": "crystal_ligand_bbox", "center_x": cx, "center_y": cy, "center_z": cz,
            "size_x": sx, "size_y": sy, "size_z": sz, "pad_A": pad_A, "max_edge_A": max_edge_A}


def prep_ligands_dekois(tg_dir: Path, cache_dir: Path, sample: dict, log=print) -> tuple[Path, dict]:
    """按采样索引把 actives/decoys 转 pdbqt；文件名 a_<idx>/d_<idx> 标注身份。"""
    out_dir = cache_dir / "lig"
    out_dir.mkdir(parents=True, exist_ok=True)
    mk = MoleculePreparation()
    stats = {"ok": 0, "skip_parse": 0, "skip_3d": 0, "skip_meeko": 0}

    def handle(kind: str, src: Path, idxs: list[int]):
        supplier = (Chem.SDMolSupplier(str(src)) if src.suffix == ".sdf"
                    else Chem.ForwardSDMolSupplier(gzip.open(src, "rb")))
        want = set(idxs)
        for i, mol in enumerate(supplier):
            if i not in want:
                continue
            dest = out_dir / f"{kind}_{i:05d}.pdbqt"
            if dest.is_file():
                stats["ok"] += 1
                continue
            if mol is None or mol.GetNumAtoms() == 0:
                stats["skip_parse"] += 1
                continue
            if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
                stats["skip_3d"] += 1
                continue
            try:
                # DEKOIS SDF 是隐式氢，Meeko 要求显式氢（PITFALLS #13）
                if any(a.GetTotalNumHs() for a in mol.GetAtoms()):
                    mol = Chem.AddHs(mol, addCoords=True)
                pdbqt, ok, _ = PDBQTWriterLegacy.write_string(mk.prepare(mol)[0])
                if not ok:
                    stats["skip_meeko"] += 1
                    continue
            except Exception:
                stats["skip_meeko"] += 1
                continue
            dest.write_text(pdbqt)
            stats["ok"] += 1

    handle("a", tg_dir / "actives.sdf", sample["actives"])
    handle("d", tg_dir / "decoys.sdfgz", sample["decoys"])
    log(f"[prep] {tg_dir.name}: {stats}")
    return out_dir, stats


def prep_ligands(set_name: str, tg: str, cache_dir: Path, sample: dict, log=print) -> tuple[Path, dict]:
    """把采样到的 actives/decoys 转成 pdbqt。文件名前缀 a_/d_ 标注身份。"""
    out_dir = cache_dir / "lig"
    out_dir.mkdir(parents=True, exist_ok=True)
    mk = MoleculePreparation()
    stats = {"ok": 0, "skip_parse": 0, "skip_3d": 0, "skip_meeko": 0}
    for kind, ids in (("a", sample["actives"]), ("d", sample["decoys"])):
        want = set(ids)
        src = registry.DATA / set_name / tg / ("actives_final.mol2.gz" if kind == "a" else "decoys_final.mol2.gz")
        for name, block in registry.read_mol2_gz(src):
            if name not in want:
                continue
            dest = out_dir / f"{kind}_{name}.pdbqt"
            if dest.is_file():
                stats["ok"] += 1
                continue
            mol = Chem.MolFromMol2Block(block, sanitize=True, removeHs=False)
            if mol is None or mol.GetNumAtoms() == 0:
                stats["skip_parse"] += 1
                continue
            if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
                stats["skip_3d"] += 1
                continue
            try:
                prepared = mk.prepare(mol)
                pdbqt, ok, _ = PDBQTWriterLegacy.write_string(prepared[0])
                if not ok:
                    stats["skip_meeko"] += 1
                    continue
            except Exception:
                stats["skip_meeko"] += 1
                continue
            dest.write_text(pdbqt)
            stats["ok"] += 1
    log(f"[prep] {tg}: {stats}")
    return out_dir, stats

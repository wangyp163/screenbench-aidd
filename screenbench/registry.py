"""冻结基准注册表：manifest + 校验和 + 采样规则。

冻结的含义：manifest 一经生成，靶点文件以 sha256 锁定，采样以 seed 锁定；
之后任何 run 的对比都基于同一套题目，分数才有可比性。
"""
import gzip
import hashlib
import json
import random
import subprocess
import time
from pathlib import Path

from .paths import REPO

DATA = REPO / "data" / "frozen"
DEKOIS_ROOT = REPO / "data" / "dekois_airdd" / "DEKOIS"
DUD_E_FILES = ["{tg}_pdb.pdb.gz", "{tg}_lig.mol2.gz", "actives_final.mol2.gz", "decoys_final.mol2.gz",
               "actives_final.ism", "decoys_final.ism"]
DEKOIS_FILES = ["ligand.sdf", "receptor.pdb", "actives.sdf", "decoys.sdfgz"]
DEFAULT_SAMPLING = {"seed": 42, "n_actives": 100, "n_decoys": 600}
DEKOIS_SAMPLING = {"seed": 42, "n_actives": 60, "n_decoys": 600}
BOX_RULE = {"rule": "crystal_ligand_bbox", "pad_A": 4.0, "max_edge_A": 26.0}


def _count_sdf(path: Path) -> int:
    """分子数 = 记录分隔符 $$$$. 行数（对 sdf 与 sdf.gz 均成立，与 RDKit supplier 索引一致）。"""
    import gzip as _gz
    op = _gz.open if path.suffix == ".sdfgz" else open
    with op(path, "rt", errors="ignore") as f:
        return sum(1 for line in f if line.rstrip("\n") == "$$$$")


def build_manifest_dekois(set_name: str = "dekois-airdd", targets: list[str] | None = None) -> Path:
    """冻结 DEKOIS(THU airdd 布局) 基准：逐文件 sha256 + 采样规则锁死。"""
    tgs = targets or sorted(p.name for p in DEKOIS_ROOT.iterdir() if p.is_dir())
    manifest = {
        "set": set_name, "layout": "dekois-airdd", "version": 2,
        "source": "hf-mirror THU-ATOM/dekois_airdd (DEKOIS 2.0, 81 targets)",
        "targets": {}, "sampling": DEKOIS_SAMPLING, "box": BOX_RULE,
    }
    for tg in tgs:
        manifest["targets"][tg] = {}
        for fname in DEKOIS_FILES:
            manifest["targets"][tg][fname] = _sha256(DEKOIS_ROOT / tg / fname)
    path = DATA / set_name / "manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    return path


def sample_indices_dekois(set_name: str, tg: str, n_actives: int | None = None,
                          n_decoys: int | None = None) -> dict:
    """按分子索引确定性采样（与 $$$$ 计数/RDKit supplier 索引一致）。

    n_* 覆盖参数用于 sweep 等场景；覆盖会被记录进报告（不静默）。
    注意：不同 k 的 sample 结果不是嵌套子集。
    """
    m = load_manifest(set_name)
    s = m["sampling"]
    rng = random.Random(s["seed"])
    want = {"actives": n_actives or s["n_actives"], "decoys": n_decoys or s["n_decoys"]}
    out = {}
    for kind, fname, key in (("a", "actives.sdf", "actives"), ("d", "decoys.sdfgz", "decoys")):
        n = _count_sdf(DEKOIS_ROOT / tg / fname)
        out[key] = sorted(rng.sample(range(n), min(want[key], n)))
    out["seed"] = s["seed"]
    out["override"] = {k: v for k, v in want.items() if v != s.get(f"n_{k}")}
    out["n_available"] = {"actives": _count_sdf(DEKOIS_ROOT / tg / "actives.sdf"),
                          "decoys": _count_sdf(DEKOIS_ROOT / tg / "decoys.sdfgz")}
    return out


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_target(set_name: str, tg: str, retries: int = 4) -> dict:
    """从 DUD-E 官方源下载靶点文件并生成带校验和的 manifest 条目。

    坑（见 PITFALLS.md #6）：DUD-E 会按 IP 限流返回 503，且 503 页面是 HTML——
    下载后必须校验内容不是 HTML，否则会把错误页冻进基准集。
    """
    tg_dir = DATA / set_name / tg
    tg_dir.mkdir(parents=True, exist_ok=True)
    entry = {}
    for tpl in DUD_E_FILES:
        fname = tpl.format(tg=tg)
        url = f"https://dude.docking.org/targets/{tg}/{fname}"
        out = tg_dir / fname
        for attempt in range(retries):
            if out.is_file() and out.stat().st_size > 1000 and not out.read_bytes().startswith(b"<"):
                break
            subprocess.run(["curl", "-sL", "-A", "Mozilla/5.0 screenbench/0.1",
                            "--max-time", "300", "-o", str(out), url])
            time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError(f"下载失败（重试 {retries} 次）: {url}")
        entry[fname] = _sha256(out)
    return entry


def build_manifest(set_name: str, targets: list[str], sampling: dict | None = None) -> Path:
    manifest = {
        "set": set_name,
        "source": "https://dude.docking.org/targets/{tg}/{file}",
        "targets": {},
        "sampling": {**DEFAULT_SAMPLING, **(sampling or {})},
        "box": BOX_RULE,
    }
    for tg in targets:
        manifest["targets"][tg] = download_target(set_name, tg)
    path = DATA / set_name / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    return path


def load_manifest(set_name: str) -> dict:
    return json.loads((DATA / set_name / "manifest.json").read_text())


def verify(set_name: str) -> list[str]:
    """校验冻结集完整性：任何文件哈希不符即报出（gate 的前置）。"""
    m = load_manifest(set_name)
    root = DEKOIS_ROOT if m.get("layout") == "dekois-airdd" else DATA / set_name
    bad = []
    for tg, files in m["targets"].items():
        for fname, want in files.items():
            p = root / tg / fname
            if not p.is_file():
                bad.append(f"{tg}/{fname}: 文件缺失")
                continue
            if _sha256(p) != want:
                bad.append(f"{tg}/{fname}: sha256 mismatch")
    return bad


def sample_ids(set_name: str, tg: str) -> dict:
    """确定性采样：同 seed 同题目。actives/decoys 各取 ism 里的 ID 列表后抽样。"""
    m = load_manifest(set_name)
    s = m["sampling"]
    rng = random.Random(s["seed"])

    def ids_of(fname):
        with open(DATA / set_name / tg / fname) as f:
            cols = [line.split() for line in f if line.strip()]
        return [c[-1] for c in cols], [c[0] for c in cols]

    act_ids, _ = ids_of("actives_final.ism")
    dec_ids, _ = ids_of("decoys_final.ism")
    rng.shuffle(act_ids)
    rng.shuffle(dec_ids)
    return {
        "actives": sorted(act_ids[: s["n_actives"]]),
        "decoys": sorted(dec_ids[: s["n_decoys"]]),
        "seed": s["seed"],
    }


def read_mol2_gz(path: Path):
    """逐条读多记录 mol2.gz，产出 (name, molblock)。name 取 @<MOLECULE> 块首行。"""
    with gzip.open(path, "rt", errors="ignore") as f:
        text = f.read()
    for block in text.split("@<TRAILER>"):
        if "@<MOLECULE>" not in block:
            continue
        mol_block = "@" + block[block.index("@<MOLECULE>"):] + "\n@<TRAILER>\n"
        lines = mol_block.splitlines()
        name = ""
        try:
            name = lines[lines.index("@<MOLECULE>") + 1].strip()
        except (ValueError, IndexError):
            pass
        yield name, mol_block

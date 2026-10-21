"""引擎适配器：统一 JSON 契约 → vina 批量 docking → 结构化结果。

契约（calcflow 味道）：入参 = 受体 + pdbqt 目录 + box + 引擎参数（全部记录）；
出参 = 每个 ligand 的 best affinity，解析自 *_out.pdbqt 的 REMARK VINA RESULT。
并行策略：把配体分片到多个 vina 进程（--cpu 1），比单进程 batch 更能吃满多核。
"""
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .paths import tool_bin, tool_env

RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)")


def _shard(paths: list[Path], n: int) -> list[list[Path]]:
    return [paths[i::n] for i in range(n)] or [[]]


def run_vina(receptor: Path, lig_dir: Path, out_dir: Path, box: dict,
             seed: int, exhaustiveness: int, n_parallel: int, log=print) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    ligs = sorted(lig_dir.glob("*.pdbqt"))
    todo = [p for p in ligs if not (out_dir / f"{p.stem}_out.pdbqt").is_file()]
    log(f"[vina] 共 {len(ligs)} 配体，待跑 {len(todo)}，并行 {n_parallel} 进程，exhaustiveness={exhaustiveness}")
    shards = _shard(todo, n_parallel)
    box_keys = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")

    def one_shard(i, shard):
        if not shard:
            return 0
        sdir = out_dir / f"shard{i}"
        sdir.mkdir(exist_ok=True)
        cmd = [tool_bin("vina"), "--receptor", str(receptor),
               "--batch", *[str(p) for p in shard], "--dir", str(sdir),
               "--cpu", "1", "--seed", str(seed), "--exhaustiveness", str(exhaustiveness)]
        for k in box_keys:
            cmd += [f"--{k}", str(box[k])]
        r = subprocess.run(cmd, capture_output=True, text=True, env=tool_env())
        if r.returncode != 0:
            log(f"[vina] shard{i} 返回码 {r.returncode}: {r.stderr[-300:]}")
        n = 0
        for f in sdir.glob("*_out.pdbqt"):
            target = out_dir / f.name
            if not target.is_file():
                target.write_bytes(f.read_bytes())
            n += 1
        return n

    with ThreadPoolExecutor(max_workers=n_parallel) as ex:
        counts = list(ex.map(lambda a: one_shard(*a), enumerate(shards)))
    return {"ligands_total": len(ligs), "docked_this_run": sum(counts)}


def parse_results(out_dir: Path, keep: set[str] | None = None) -> dict[str, float]:
    """解析各配体 best-mode affinity。keep：只保留指定 ligand 名（采样隔离，防缓存串味）。"""
    scores = {}
    for f in sorted(out_dir.glob("*_out.pdbqt")):
        stem = f.stem[:-4] if f.stem.endswith("_out") else f.stem
        if keep is not None and stem not in keep:
            continue
        m = RESULT_RE.search(f.read_text(errors="ignore"))
        if m:
            scores[stem] = float(m.group(1))
    return scores


def engine_version() -> str:
    r = subprocess.run([tool_bin("vina"), "--version"], capture_output=True, text=True, env=tool_env())
    return (r.stdout or r.stderr).strip().splitlines()[0]

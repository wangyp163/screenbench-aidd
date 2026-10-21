"""screenbench CLI：manifest-build / verify / run / gate / audit"""
import argparse
import time
from pathlib import Path

from . import registry, prep, engines, metrics, report, audit
from .paths import REPO

RUNS = REPO / "runs"


def cmd_manifest(a):
    if a.layout == "dekois-airdd":
        p = registry.build_manifest_dekois(a.set, a.targets)
    else:
        p = registry.build_manifest(a.set, a.targets)
    print(f"manifest 写入 {p}")


def _run_dekois(a, m):
    run_dir = RUNS / f"{a.set}_{a.engine}_{time.strftime('%Y%m%d_%H%M%S')}"
    run_dir.mkdir(parents=True, exist_ok=True)
    targets_out = {}
    for tg in (a.targets or list(m["targets"])):
        tg_dir = registry.DEKOIS_ROOT / tg
        cache = REPO / "cache" / a.set / tg
        t0 = time.time()
        rec = prep.prep_receptor_dekois(tg_dir, cache)
        box = prep.box_from_ligand_sdf(tg_dir, **{k: v for k, v in m["box"].items() if k != "rule"})
        sample = registry.sample_indices_dekois(a.set, tg)
        lig_dir, stats = prep.prep_ligands_dekois(tg_dir, cache, sample)
        out_dir = cache / "dock" / f"exh{a.exhaustiveness}"
        run = engines.run_vina(rec, lig_dir, out_dir, box, a.seed, a.exhaustiveness, a.parallel)
        scores = engines.parse_results(out_dir)
        labels = [1 if k.startswith("a_") else 0 for k in scores]
        targets_out[tg] = metrics.evaluate(labels, list(scores.values()))
        skipped = stats.get("skip_parse", 0) + stats.get("skip_3d", 0) + stats.get("skip_meeko", 0)
        print(f"[run] {tg}: {targets_out[tg]}  用时 {time.time()-t0:.0f}s  skip={skipped}")
    payload = {
        "set": a.set, "engine": a.engine, "engine_version": engines.engine_version(),
        "engine_params": {"seed": a.seed, "exhaustiveness": a.exhaustiveness, "parallel": a.parallel},
        "sampling": {"seed": sample["seed"], "n_available": sample["n_available"],
                     "n_sampled": {"actives": len(sample["actives"]), "decoys": len(sample["decoys"])}},
        "box": box, "targets": targets_out, "prep_stats": stats, "docked": run,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    return run_dir, payload


def cmd_verify(a):
    bad = registry.verify(a.set)
    print("冻结集校验通过" if not bad else "校验失败:\n" + "\n".join(bad))
    raise SystemExit(0 if not bad else 1)


def cmd_run(a):
    m = registry.load_manifest(a.set)
    if m.get("layout") == "dekois-airdd":
        run_dir, payload = _run_dekois(a, m)
    else:
        run_dir, payload = _run_dude(a, m)
    j, md = report.write_report(run_dir, payload)
    print(f"报告: {md}")
    if a.gate:
        ok, failures = report.check_gate(j, a.gate)
        print("gate: PASS" if ok else "gate: FAIL\n" + "\n".join(failures))
        raise SystemExit(0 if ok else 1)


def _run_dude(a, m):
    run_dir = RUNS / f"{a.set}_{a.engine}_{time.strftime('%Y%m%d_%H%M%S')}"
    run_dir.mkdir(parents=True, exist_ok=True)
    targets_out = {}
    for tg in (a.targets or list(m["targets"])):
        cache = REPO / "cache" / a.set / tg
        t0 = time.time()
        rec = prep.prep_receptor(a.set, tg, cache)
        box = prep.box_from_crystal_ligand(a.set, tg, **{k: v for k, v in m["box"].items() if k != "rule"})
        sample = registry.sample_ids(a.set, tg)
        lig_dir, stats = prep.prep_ligands(a.set, tg, cache, sample)
        out_dir = cache / "dock" / f"exh{a.exhaustiveness}"
        run = engines.run_vina(rec, lig_dir, out_dir, box, a.seed, a.exhaustiveness, a.parallel)
        scores = engines.parse_results(out_dir)
        labels = [1 if k.startswith("a_") else 0 for k in scores]
        targets_out[tg] = metrics.evaluate(labels, list(scores.values()))
        skipped = stats.get("skip_parse", 0) + stats.get("skip_3d", 0) + stats.get("skip_meeko", 0)
        print(f"[run] {tg}: {targets_out[tg]}  用时 {time.time()-t0:.0f}s  skip={skipped}")
    payload = {
        "set": a.set, "engine": a.engine, "engine_version": engines.engine_version(),
        "engine_params": {"seed": a.seed, "exhaustiveness": a.exhaustiveness, "parallel": a.parallel},
        "sampling": m["sampling"], "box": box, "targets": targets_out,
        "prep_stats": stats, "docked": run,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    return run_dir, payload


def cmd_audit(a):
    print(audit.audit_csv(a.csv, a.col))


def main():
    p = argparse.ArgumentParser(prog="screenbench")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("manifest-build"); s.add_argument("--set", required=True)
    s.add_argument("--layout", default="dekois-airdd", choices=["dekois-airdd", "dude"])
    s.add_argument("--targets", nargs="+"); s.set_defaults(fn=cmd_manifest)

    s = sub.add_parser("verify"); s.add_argument("--set", required=True); s.set_defaults(fn=cmd_verify)

    s = sub.add_parser("run"); s.add_argument("--set", required=True)
    s.add_argument("--engine", default="vina", choices=["vina"])
    s.add_argument("--targets", nargs="+")
    s.add_argument("--seed", type=int, default=42)
    s.add_argument("--exhaustiveness", type=int, default=4)
    s.add_argument("--parallel", type=int, default=13)
    s.add_argument("--gate", nargs="*"); s.set_defaults(fn=cmd_run)

    s = sub.add_parser("audit"); s.add_argument("--csv", required=True)
    s.add_argument("--col", default="smiles"); s.set_defaults(fn=cmd_audit)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()

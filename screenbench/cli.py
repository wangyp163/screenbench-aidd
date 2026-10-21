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
    targets_out, boxes, last_stats = {}, {}, {}
    for tg in (a.targets or list(m["targets"])):
        tg_dir = registry.DEKOIS_ROOT / tg
        cache = REPO / "cache" / a.set / tg
        t0 = time.time()
        rec = prep.prep_receptor_dekois(tg_dir, cache)
        box = prep.box_from_ligand_sdf(tg_dir, **{k: v for k, v in m["box"].items() if k != "rule"})
        boxes[tg] = box
        sample = registry.sample_indices_dekois(a.set, tg, a.n_actives, a.n_decoys)
        lig_dir, stats = prep.prep_ligands_dekois(tg_dir, cache, sample)
        last_stats = stats
        # 采样隔离缓存：不同采样配置的结果不混目录（否则 parse 会串味，PITFALLS #17）
        tag = f"exh{a.exhaustiveness}_a{len(sample['actives'])}_d{len(sample['decoys'])}"
        out_dir = cache / "dock" / tag
        run = engines.run_vina(rec, lig_dir, out_dir, box, a.seed, a.exhaustiveness, a.parallel)
        keep = {f"a_{i:05d}" for i in sample["actives"]} | {f"d_{i:05d}" for i in sample["decoys"]}
        scores = engines.parse_results(out_dir, keep=keep)
        labels = [1 if k.startswith("a_") else 0 for k in scores]
        targets_out[tg] = metrics.evaluate(labels, list(scores.values()))
        skipped = stats.get("skip_parse", 0) + stats.get("skip_3d", 0) + stats.get("skip_meeko", 0)
        print(f"[run] {tg}: {targets_out[tg]}  用时 {time.time()-t0:.0f}s  skip={skipped}", flush=True)
    payload = {
        "set": a.set, "engine": a.engine, "engine_version": engines.engine_version(),
        "engine_params": {"seed": a.seed, "exhaustiveness": a.exhaustiveness, "parallel": a.parallel},
        "sampling": {"seed": sample["seed"], "n_available": sample["n_available"],
                     "n_sampled": {"actives": len(sample["actives"]), "decoys": len(sample["decoys"])},
                     **({"override": sample["override"]} if sample.get("override") else {})},
        "box": boxes, "targets": targets_out, "prep_stats": last_stats, "docked": run,
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


def cmd_gate(a):
    ok, failures = report.check_gate(Path(a.report), a.rules)
    print("gate: PASS" if ok else "gate: FAIL\n" + "\n".join(failures))
    raise SystemExit(0 if ok else 1)


def cmd_redock(a):
    from . import redock
    redock.redock_and_report(a.set, a.target, a.exhaustiveness, a.threshold)


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
    s.add_argument("--n-actives", type=int, default=None, help="覆盖采样(记录进报告)")
    s.add_argument("--n-decoys", type=int, default=None, help="覆盖采样(记录进报告)")
    s.add_argument("--gate", nargs="*"); s.set_defaults(fn=cmd_run)

    s = sub.add_parser("gate"); s.add_argument("--report", required=True)
    s.add_argument("--assert", dest="rules", nargs="+", required=True,
                   metavar="RULE", help='如 "all:ROC_AUC>=0.5" "80%:EF5_pct>=2"')
    s.set_defaults(fn=cmd_gate)

    s = sub.add_parser("redock"); s.add_argument("--set", default="dekois-airdd")
    s.add_argument("--target", required=True)
    s.add_argument("--exhaustiveness", type=int, default=8)
    s.add_argument("--threshold", type=float, default=2.0)
    s.set_defaults(fn=cmd_redock)

    s = sub.add_parser("audit"); s.add_argument("--csv", required=True)
    s.add_argument("--col", default="smiles"); s.set_defaults(fn=cmd_audit)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()

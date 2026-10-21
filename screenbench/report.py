"""报告与门禁：JSON 契约输出 + markdown 摘要 + 断言 gate。

gate 规则形如 "EF5_pct>=3"（对全部靶点生效）或 "80%:ROC_AUC>=0.55"（80% 靶点达标即可）。
机检不过 = run 不可用于结论，防止挑靶点报数（这正是本工具存在的理由）。
"""
import json
import platform
from datetime import datetime
from pathlib import Path

from .paths import REPO


def write_report(run_dir: Path, payload: dict) -> tuple[Path, Path]:
    j = run_dir / "report.json"
    j.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    lines = [f"# screenbench 报告 · {payload['set']} · {payload['engine']}",
             f"- 时间：{payload['timestamp']}　引擎：{payload['engine_version']}　host: {platform.node()}",
             f"- 采样：{payload['sampling']}　box：{payload['box']}",
             f"- exhaustiveness={payload['engine_params']['exhaustiveness']} seed={payload['engine_params']['seed']}",
             "", "| 靶点 | n(active/total) | EF1% | EF5% | ROC-AUC |", "|---|---|---|---|---|"]
    for tg, m in payload["targets"].items():
        lines.append(f"| {tg} | {m['n_active']}/{m['n_total']} | {m['EF1_pct']} | {m['EF5_pct']} | {m['ROC_AUC']} |")
    prep = payload.get("prep_stats", {})
    if prep:
        lines += ["", f"前处理统计：{prep}"]
    lines += ["", "> 冻结集 + 固定 seed，任何复跑应逐字段一致（除时间戳/主机名）。", ""]
    md = run_dir / "report.md"
    md.write_text("\n".join(lines))
    return j, md


def check_gate(report_json: Path, asserts: list[str]) -> tuple[bool, list[str]]:
    payload = json.loads(report_json.read_text())
    failures = []
    for rule in asserts:
        quant, expr = (rule.split(":", 1) + [None])[:2] if ":" in rule else ("all", rule)
        if expr is None:
            quant, expr = "all", rule
        metric, _, cond = expr.partition(">=")
        threshold = float(cond)
        values = {tg: m[metric] for tg, m in payload["targets"].items()}
        passed = {tg: v for tg, v in values.items() if isinstance(v, (int, float)) and v >= threshold}
        if quant == "all":
            ok = len(passed) == len(values)
        elif quant.endswith("%"):
            ok = len(passed) / max(1, len(values)) >= int(quant[:-1]) / 100
        else:
            ok = False
        if not ok:
            failures.append(f"{rule} → 达标 {len(passed)}/{len(values)} {values}")
    return (not failures), failures

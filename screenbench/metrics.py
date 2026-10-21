"""指标：EF1%/EF5%（富集因子）与 ROC-AUC。

EF 的定义：top-K 里 active 占比 / 全库 active 占比，K = round(fraction * N)。
注意：子采样集上 EF1% 统计功效弱（见 PITFALLS.md），因此 EF5% 一并报告。
"""
from sklearn.metrics import roc_auc_score


def enrichment_factor(labels: list[int], scores: list[float], fraction: float) -> float:
    n = len(labels)
    k = max(1, round(fraction * n))
    order = sorted(range(n), key=lambda i: scores[i], reverse=True)
    top = order[:k]
    n_act = sum(labels)
    if n_act == 0:
        return float("nan")
    hit_rate_top = sum(labels[i] for i in top) / k
    base_rate = n_act / n
    return hit_rate_top / base_rate


def evaluate(labels: list[int], scores: list[float], score_dir: str = "lower_better") -> dict:
    """score_dir: docking 打分约定越小越好（vina 负值），内部统一翻转成越大越好再算指标。"""
    s = [-v for v in scores] if score_dir == "lower_better" else list(scores)
    return {
        "EF1_pct": round(enrichment_factor(labels, s, 0.01), 3),
        "EF5_pct": round(enrichment_factor(labels, s, 0.05), 3),
        "ROC_AUC": round(roc_auc_score(labels, s), 4),
        "n_total": len(labels),
        "n_active": int(sum(labels)),
    }

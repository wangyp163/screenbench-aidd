# screenbench 报告 · dekois-airdd · vina
- 时间：2026-10-21 (metrics recomputed from cached poses; sign bug fixed, PITFALLS #16)　引擎：AutoDock Vina v1.2.5　host: wangyp-pc
- 采样：{'seed': 42, 'n_sampled': {'actives': 60, 'decoys': 600}}　box：{'rule': 'crystal_ligand_bbox pad=4A max=26A'}
- exhaustiveness=4 seed=42

| 靶点 | n(active/total) | EF1% | EF5% | ROC-AUC |
|---|---|---|---|---|
| cox2 | 60/660 | 9.429 | 10.333 | 0.9539 |
| ache | 60/660 | 4.714 | 3.333 | 0.6792 |
| adrb2 | 60/660 | 3.143 | 3.333 | 0.6565 |

前处理统计：{'ok': 660, 'skip': 0}

> 冻结集 + 固定 seed，任何复跑应逐字段一致（除时间戳/主机名）。

# screenbench

**虚拟筛选统一评测 harness** —— 冻结基准 + 统一指标 + 断言门禁 + 数据审计。
一句话：*pytest for virtual screening*。

## 为什么需要它

虚拟筛选领域的日常：各自挑靶点、各自调盒子、各自报好看的数字，跨论文无法比较。
screenbench 把评测变成工程契约：

- **冻结基准**：靶点数据以 sha256 锁进 manifest，采样以 seed 锁死——同一把尺子量所有引擎
- **统一指标**：EF1% / EF5% / ROC-AUC（显式声明打分方向，见 PITFALLS #16）
- **断言门禁**：`gate` 规则不通过 = 结论不可用，从机制上堵死 cherry-picking
- **数据审计**：任意分子 CSV 的骨架多样性 + 集内近重复检测

## 首批数字（v0.1 · 2026-10-20）

AutoDock Vina 1.2.5 · DEKOIS 2.0 · 每靶点采样 60 actives + 600 decoys（seed=42）· exh=4 · 16 核笔记本

| 靶点 | EF1% | EF5% | ROC-AUC |
|---|---|---|---|
| cox2 | 9.43 | 10.33 | **0.954** |
| ache | 4.71 | 3.33 | 0.679 |
| adrb2 | 3.14 | 3.33 | 0.657 |

> 冻结集 + 固定 seed，复跑应逐字段一致（报告在 `runs/`）。cox2 是 Vina 公认好打的靶点，ache/adrb2 偏难——方向与文献一致。

## 快速开始

```bash
# 依赖：vina + openbabel（tools/ 下 deb 本地解包，免 sudo，见 PITFALLS.md）
pip install -e .
screenbench manifest-build --set dekois-airdd --targets cox2 ache adrb2  # 冻结基准
screenbench verify --set dekois-airdd                                    # sha256 校验
screenbench run --set dekois-airdd --engine vina --seed 42 --exhaustiveness 4
screenbench audit --csv your_library.csv                                 # 泄漏/近重复审计
```

## 组件

| 组件 | 职责 | 状态 |
|---|---|---|
| registry | 冻结 manifest（sha256 + 采样 seed + 盒子规则） | ✅ |
| engines | Vina 适配器（声明式契约，分片并行批处理） | ✅ vina |
| metrics | EF1%/EF5%/AUC，显式打分方向 | ✅ |
| report/gate | JSON+MD 报告、断言门禁 | ✅ |
| audit | 骨架多样性 + Tanimoto 近重复 | ✅ |
| redocking RMSD 验证 | 晶体配体重对接 <2Å 证明管线正确性 | 🔜 v0.2 |
| 更多引擎 / 全 81 靶点 / DUD-E 布局 | smina、AutoDock-GPU | 🔜 |

## 开发实录

- **[PITFALLS.md](PITFALLS.md)**：16 个真实踩坑（DUD-E 503 返回 HTML 错误页、hf-mirror 限速 40 倍绕行、打分方向反了的首轮翻车……）——评测基建为什么难，这里全是证据
- **[PROJECT.md](PROJECT.md)**：设计蓝图与定位
- `phase0_split_leakage.py`：项目的前身实验——random vs scaffold split 泄漏实测（随机切分虚高 2.6%，为什么冻结评测集必须配防泄漏切分）

## License

MIT

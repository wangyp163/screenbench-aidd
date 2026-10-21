# screenbench · 项目蓝图（v0 草案，2026-10-20）

> 定位：虚拟筛选（virtual screening）的统一评测 harness——冻结基准 + 统一指标 + 断言门禁 + 泄漏审计。
> 一句话：**pytest for virtual screening**。
> 为什么是它：CADD 开源圈模型的人多、评测的人几乎没有；冻结评测集/断言/回归门禁这套方法论是我的差异化，且直接命中 AIDD 岗 C（结构辅助评价流程）+ E（统一评测、研究组件复用）。

## CLI 形态

```
screenbench registry list                                # 冻结基准注册表（带校验和/manifest）
screenbench run --set dude-v1 --engine vina              # 引擎适配器
screenbench report runs/<id>.json                        # EF1%/EF5%、ROC-AUC、BedROC、RMSD 成功率
screenbench gate runs/<id>.json --assert "EF1%>=20 on 80% targets"   # CI 门禁
screenbench audit --dataset my.csv                       # 泄漏审计：scaffold 分布 + Tanimoto 近重复
```

## 组件 × 已有资产映射

| 组件 | v0.1 内容 | 复用资产 |
|---|---|---|
| 冻结基准注册表 | DUD-E 3 靶点子集，manifest+校验和，版本化 | 压测计划"先确认再执行" |
| 引擎适配器 | Vina（CPU 并行）；声明式 JSON 契约 | calcflow pre/run/post 协议 |
| 指标层 | EF1%/EF5%、ROC-AUC、redocking RMSD(<2Å) 成功率 | DUD-E actives/decoys 金标 |
| 断言门禁 | 配置驱动阈值、与基线回归对比 | 压测 SOP + stress-report-lint |
| 泄漏审计 | scaffold 切分 + 训练/测试近重复检测 | Phase 0 实验直接移植 |

## v0.1 排期（两个周末）

- 周末 1：Vina 适配器 + 3 靶点跑通 + EF/AUC + JSON 报告 + manifest
- 周末 2：gate 断言 + audit 子命令 + README（实测数字开场）+ GitHub Actions mini 回归
- v0.2+：Smina/AutoDock-GPU 适配器、LIT-PCBA（低偏差靶点集）、BedROC、失败案例库（高分 decoy 归档）

## 运营

- README 第一屏放 Vina×3 靶点 EF1% 实测表
- MIT；贡献点友好（加一个靶点=一个 PR，terminal-bench 众包模式同构）
- 冷启动：r/comp_chem 帖（主题：冻结基准下 cherry-picking 的代价，用 2.6% 泄漏数字开场）+ RDKit mailing list + 智药邦

## 立场与红线

- 本项目是**评测基建**，不训模型、不吹 SOTA；所有数字可复现（冻结集 + 固定种子 + 版本钉扎）
- 与 demo-ai4s（材料域）并列，构成"同一套评测方法论跨域"的作品集叙事

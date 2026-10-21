# PITFALLS · 坑与优化点记录（screenbench 开发全程）

> 目的：真实踩坑流水账——既是项目复盘素材，也是"评测基建为什么难"的实证。按时间顺序追加。

## 环境与安装

**#1 PyPI 直连极慢，装个 rdkit 280s 超时**
- 现象：`pip install rdkit` 直连 4 分钟无输出，pip 进程消失后 site-packages 为空
- 解法：清华镜像 `-i https://pypi.tuna.tsinghua.edu.cn/simple`，秒级
- 教训：CI/文档里必须写默认镜像逻辑或 `PIP_INDEX_URL` 检查

**#2 PyTDC 不兼容 Python 3.13（metadata-generation-failed）**
- 现象：`pip install PyTDC` 在 py3.13 报子进程元数据错误（老包 setup.py 链路）
- 解法：不装 PyTDC，数据集从源头直下（MoleculeNet ↔ deepchemdata S3）
- 优化点：screenbench 的 `audit` 子命令本身就该取代 PyTDC 的数据获取角色

**#3 pip 版 vina 没有 cp313 wheel**
- 现象：`pip install vina` 失败（官方 wheel 只到 3.12）
- 解法：apt 下载 deb + `dpkg -x` 本地解包（免 sudo）：
  `apt-get download autodock-vina openbabel libopenbabel7 && dpkg -x *.deb tools/local`
- 附带坑 #3b：vina 二进制缺 `libboost_{filesystem,program_options}.so.1.83.0` → 再下载两个 boost deb 解包，运行时注入 `LD_LIBRARY_PATH`（见 `screenbench/paths.py: tool_env()`）
- 设计决策：引擎走"外部二进制 + CLI 契约"而不是 Python API——正好符合"引擎适配器"架构，还绕开了版本地狱

## 数据获取

**#4 DUD-E 官方源 503（最大的坑）**
- 现象：`dude.docking.org/targets/ache/*` 四次重试全部 503；更阴险的是第一次 actives_final.ism 拿到 299 字节的 HTML 错误页（若不校验会冻进基准集）
- 分析：503 是服务端按 IP 限流/过载（历史上常年如此），换出口 IP（Clash 代理）或镜像才能绕
- 解法（三重）：① 下载后校验 magic bytes 非 `<`（HTML）+ 最小体积；② UA + 指数退避；③ 镜像兜底
- 转机：Zenodo 直连可达，**DEKOIS 2.0 pdbqt 版**（zenodo record 17050764，128.5MB 整包）——科学上还更优（DEKOIS 就是针对 DUD-E 类似物偏差设计的），且 pdbqt 现成、省掉配体制备链路
- 教训：**冻结基准的第一课就是"数据源不可靠"**——manifest 必须带 sha256 + 内容校验，这本身就是产品特性

**#5 raw.githubusercontent.com 被墙（000）但 codeload.github.com 通**
- GitHub 单文件路径要走 codeload 整包或 jsdelivr CDN；Clash 开启后此限制解除

**#6 MoleculeNet HIV.csv 下载疑似截断**
- 现象：S3 下载 1.5MB/28.6k 行（官方全量 41k），curl 无报错
- 处理：实验有效性不受影响（两种切分同一数据对比），报告按实际行数写；大文件下载必须加行数/sha256 断言

## 设计决策（预判的坑，提前规避）

**#7 vina 1.2 batch 模式的并行姿势**
- vina 单进程内多线程吃不满多核（大量小配体场景）→ 分片成 N 个进程 × `--cpu 1` 并行（默认 13 片留 1 核给系统）
- 优化点：后续按配体大小动态分片（大配体慢，会拖尾）

**#8 子采样集上 EF1% 统计功效弱**
- 600 decoys 的 top1% 只有 6 个分子，EF1% 方差巨大 → 指标层同时报 EF5%，报告里标注采样规模
- 这是"诚实评测"的一部分：指标必须和样本量匹配，否则就是新的 cherry-picking

**#9 配体 3D 用数据源自带构象（DUD-E mol2 是 CORINA 生成的），不现场 gen3d**
- 理由：确定性（同 seed 同构象）+ 速度 + 与文献可比
- 若换 SMILES-only 数据源，需要补 ETKDG+质子化，并固定随机种子（优化点：接受外部 conf 生成器注入）

**#10 盒子规则必须进 run config**
- 盒子定义（晶体配体 bbox + pad + 封顶）不同实验室差异巨大，是"隐性调参黑洞"——screenbench 把盒子规则作为冻结 manifest 的一部分，报告逐字段记录

## 性能账（首轮实测，2026-10-20，16 核 / exh=4 / seed=42）

| 靶点 | 配体数 | wall time | 吞吐 | 备注 |
|---|---|---|---|---|
| cox2 | 660 | ~13 min | ~51/min | 中等配体 |
| ache | 660 | 明显更慢（~4-18/min 波动） | — | 配体更大更柔，**印证 #7 的"按配体大小动态分片"优化点**：等分片会让大配体拖尾 |
| meeko 制备 | 660 | ~1.5 min | ~440/min | 含 AddHs |

- 冒烟（exh=1）：4 配体 2.7s，活性/非活性打分方向正确
- hf-mirror 下载：单连接 ~160KB/s → 12 连接分片 ~6.5MB/s

---

## 追记（22:00–23:00 施工过程）

**#11 pkill 匹配到自己的命令行，把自己杀了**
- 现象：`pkill -f "dekois_airdd.zip"` 时，当前 shell 的命令行里也含这个串，直接自杀，命令无输出
- 解法：先 `pgrep -af` 列出确认 PID 再精确 kill；或 pkill pattern 排除自身

**#12 obabel deb 解包版找不到插件库**
- 现象：`obabel -O out.pdbqt` 报 "Unable to find OpenBabel plugins. Try setting the BABEL_LIBDIR"
- 解法：`BABEL_LIBDIR=<tools>/local/usr/lib/x86_64-linux-gnu/openbabel/3.1.1`（已固化进 `paths.tool_env()`）

**#13 Meeko 拒绝隐式氢分子**
- 现象：`MoleculePreparation.prepare(mol)` 抛 "RDKit molecule has implicit Hs. Need explicit Hs."
- 解法：`Chem.AddHs(mol, addCoords=True)`（带坐标加氢，保持 3D 构象）
- 教训：skip_meeko 计数从 4→0 就是这么修的；数据源格式假设（sdf 一定带显式氢？不一定）不能有

**#14 hf-mirror 单连接限速 ~160KB/s，711MB 要 75 分钟**
- 解法：range 分片 + 12 并行连接，聚合 ~6.5MB/s（40 倍加速），2 分钟拉完
- 附带坑：分片中会有个别慢连接爬行（60KB/s）——杀掉重开新连接就恢复；断点续传逻辑进 fetch 脚本
- 佐证：THU-ATOM/dekois_airdd.zip（711,510,511 字节）拼装后 unzip -t 全绿

**#15 数据源阴差阳错成为更优解**
- 原计划 DUD-E（mol2，要做格式转换）→ 503 后改道 Zenodo DEKOIS pdbqt → 又发现 hf-mirror 上的 THU-ATOM/dekois_airdd（receptor.pdb + ligand.sdf + actives.sdf + decoys.sdfgz，81 靶点全量）
- 收获：后者走自家制备链路（obabel+Meeko），被测对象完整；且 DEKOIS 2.0 在科学上比 DUD-E 更严谨（无类似物偏差）
- 教训：评测工具的"基准集适配层"必须抽象（registry 已支持 layout: dekois-airdd / dude），数据源永远会出幺蛾子

**冒烟测试结果（cox2，exh=1）**：4/4 配体全链路通过，2.7s；活性分子打分 -11.1/-11.3 vs decoy -8.2/-8.1，方向正确 ✓

**#16 打分方向反了——首轮三个靶点 AUC 全部 <0.5（最贵的一课）**
- 现象：首次正式跑分 cox2 AUC=0.046、ache=0.32、adrb2=0.34，EF 全 0——"比随机还差"必然有 bug
- 根因：Vina 打分是**负值、越小越好**，而 sklearn `roc_auc_score` 和 EF 的 top-K 排序都默认**越大越属于正类**——方向整体反了
- 讽刺的是冒烟数据早把答案写在墙上：active -11.1 vs decoy -8.2，"越负越好"一目了然
- 修复：`metrics.evaluate` 增加 `score_dir="lower_better"` 参数，内部统一翻转（cox2 修正后 AUC=0.954）
- 教训：**指标层必须显式声明打分方向**；"结果差得离谱"永远是查 bug 的信号，不是调参的信号。此坑直接进 README 的设计原则

---

## 第二日（2026-10-21）：v0.2 施工

**#17 采样变更会污染缓存（设计层坑）**
- 场景：sweep 用 30/300 采样，但 lig pdbqt 缓存和 dock 输出目录与 60/600 共用——`parse_results` 会把旧采样的大分子集一起读进来，报告声称 30/300 实际算的是 660 个
- 修复：① dock 输出目录带采样标签（`exh4_a30_d300`）；② `parse_results` 增加 `keep` 参数按采样名单过滤
- 教训：**缓存的键必须包含所有影响语义的参数**，否则"缓存加速"变成"缓存串味"

**#18 DUD-E 第二基准的三条路一晚全灭**
- 官方直链：从 503 变 404——服务器端文件布局已改，下载走 /generate 表单（不可脚本化）
- Figshare 镜像（4.32GB）：API 和页面都 403（挡数据中心/代理 IP）
- jsdelivr 镜像 biobricks-ai/dude：文件列表 API 超时（仓库太大）
- 结论：DUD-E 收编延后；DEKOIS(81 靶点) 为主基准足够支撑 v0.x 叙事
- 优化点：找一个可脚本化的 DUD-E 子集源（社区 mirror 待调研）

**#19 verify 的路径 bug + shell 的 `&&`+`&` 组合拳**
- `verify` 沿用 DUD-E 的 `data/frozen/<set>/<tg>` 路径，但 dekois 布局实际在 `data/dekois_airdd/DEKOIS/<tg>`——manifest 建得了、校验找不到文件
- 更蠢的是启动命令 `verify && nohup ... & echo`：`&` 把整个 `&&` 链丢后台，verify 失败短路后 nohup 根本没跑，还打印了假的 PID
- 教训：后台启动要单独成行；"启动成功"必须以进程/产物为证，不能信 echo

**redocking 验证（v0.2 收官件）✅**
- cox2：RMSD **1.12 Å**（阈值 2 Å），affinity -11.25，exh=8 单配体 15.5s
- 整条制备+docking 管线（obabel 受体 + Meeko 配体 + 盒子规则 + vina）自此有正确性证明
- 实现注意：RMSD 不绕 sdf 转换——vina 输出 pdbqt 保持输入原子顺序，直接对 MODEL 1 重原子坐标逐对计算（避开 obabel 往返的拓扑损伤风险）

**全 81 靶点 sweep 已挂**：30 actives + 300 decoys/靶点（覆盖参数已显式记录进报告，不静默）、exh4、13 进程。预计 10+ 小时，缓存可断点续跑。


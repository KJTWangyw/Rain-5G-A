# Rain-5G-A

## 短距毫米波链路测雨误差预算

本分支 `paper/openmrg-error-budget-zh` 基于提交 `3174afd04993fc829858cf5fca0228dd6a351ed1`，对应论文《基于跨链路回归与量化重放的短距毫米波链路测雨误差预算》的中文科学核验修订。

研究以OpenMRG商业微波回传链路为实测对象，将湿期统计截距、晴空参考漂移、功率量化重放和理想NR导频参考噪声放在衰减域中比较。中文稿见 [paper/Manuscript_zh.pdf](paper/Manuscript_zh.pdf)，可编辑版本见 [paper/Manuscript_zh.docx](paper/Manuscript_zh.docx)。本轮PDF为核验修订稿，量化历史统计的重算状态在稿内和本页单独标记。

## 当前结果及适用范围

- 24–40 GHz筛选范围内共有660条子链路，其中646条进入分析，对应323条物理路径。实际频率为28.2–38.7 GHz，长度为0.14–11.1 km。
- 全季湿期回归截距约为0.54 dB，按物理路径簇自助得到的95%区间为0.51–0.58 dB。该量是既定基线和雷达支持条件下的统计截距，其物理组成还需要独立观测约束。
- 采用6–7月标定、8月同一网络留出，并由雷达辅助干湿分类，长度不超过0.5 km组的链路级小时雨强RMSE中位数从5.18降至1.53 mm/h，约降低70%。校正后的链路与雷达累计雨量比中位数仍为3.1，说明剩余系统高估仍需控制。
- 在792个有效导频、15 dB逐导频信噪比和10 ms相干参考的理想模型中，共享参考热噪声标准差约为0.0043 dB。真实NR接收机中的增益、波束、信道变化和幅度标定需要实测评价。

8月评价沿用雷达定义的干湿窗口，雷达确认的干期输出置零，湿期窗口包含前后30分钟的信息。因此，上述改进对应雷达辅助的同网时间留出评价。自主干湿识别、因果实时处理、独立降雨参考和跨网络迁移属于后续验证范围。

## 本轮核验与量化修复状态

| 项目 | 本轮状态 |
|---|---|
| 归档中间表复算 | `audit_summary.py`完成99项检查，当前归档核验结果为全部通过，输出见`results/budget/audit_summary.json` |
| 复算内容 | 链路筛选、物理路径簇自助区间、回归截距及斜率、分组汇总和主要留月评价指标 |
| 量化共同有效支持 | 已修复TSL、RSL及量化版本的有效样本对齐，加入3项针对分析主入口的回归测试 |
| 量化历史统计 | 参考端0.19 dB、雨期端0.15 dB、事件级合成0.14 dB及图3右图仍来自修复前归档，待原始时序重跑后核定 |
| 原始数据处理 | 原分钟缓存目前不可用，本轮原始OpenMRG至分钟缓存的处理未重跑 |
| NR时域仿真 | 本轮复核归档输出与理论口径，NR波形生成及估计样本未重新生成 |
| 中文文稿 | 已修订证据范围、统计定义、已有研究定位和排版，当前PDF标记为核验修订稿 |

归档审计从已有中间表重新计算统计量。其输入指纹、软件版本、检查结果和覆盖范围均写入`audit_summary.json`。逐事件量化RMS、原始分钟聚合、雷达路径采样、逐事件参考窗及NR波形仿真仍需通过完整流程核验。

量化修复后的结果需要恢复原始OpenMRG资料，重新执行`prepare_openmrg.py`、`analyze_openmrg.py`、`summarize.py`和`make_figures.py`，随后再次汇总并构建文稿。图3右图、量化数值及相关比较应随本次重算一并更新。

## 仓库结构

| 路径 | 内容 |
|---|---|
| `src/cmlbudget/openmrg.py` | OpenMRG读取、分钟电平、量化报告及雷达路径平均 |
| `src/cmlbudget/retrieval.py` | 干湿窗口、事件参考、湿期衰减和雨强换算 |
| `src/cmlbudget/quantization.py` | 功率量化器、高斯电平模型下的平均量化偏差及其RMS |
| `src/cmlbudget/offset.py` | 跨链路截距回归、雨强相关偏移函数和数值反演 |
| `src/nrpilot/` | NR DM-RS链路、剖面估计器、信息下界、误差传播及ITU-R雨衰关系 |
| `experiments/budget/` | 数据下载、准备、分析、导频仿真、汇总、作图及归档审计入口 |
| `experiments/budget/audit_summary.py` | 从归档中间表复算主要统计量和物理路径簇自助区间 |
| `results/budget/` | 归档分析表、`summary.json`及本轮`audit_summary.json` |
| `paper/` | 中文稿源文件、参考模板、图表、Word构建器及中文Word和PDF |
| `tests/test_quantization_support.py` | 缺失TSL、量化版本共同支持及完整功率支持的3项回归测试 |

原导频覆盖率研究的入口位于`experiments/run_main.py`、`summarize.py`、`benchmark_timing.py`和`make_figures.py`，对应结果位于`results/*.csv`。本稿误差预算流程的入口集中在`experiments/budget/`。

## 环境与核验命令

从仓库根目录执行命令。Python依赖由`requirements.txt`安装，其中`pandas>=2.2`支持汇总脚本使用的`groupby.apply(..., include_groups=False)`参数。本轮归档审计使用Python 3.12.14、NumPy 2.3.5和pandas 2.2.3，具体版本同时保存在审计输出中。

```bash
pip install -r requirements.txt
python experiments/budget/audit_summary.py
python -m unittest discover -s tests -p test_quantization_support.py -v
```

当前归档对应99项审计检查。3项回归测试使用构造数据运行分析主入口，覆盖共同有效支持修复。审计发现数值不一致时以非零退出码结束，便于接入持续集成。

## 完整复现流程

下列命令保留从下载数据到生成中文稿的完整顺序。首次运行需要下载约318 MB的OpenMRG资料，分钟缓存约占2 GB。处理耗时随硬件和软件环境变化。

```bash
pip install -r requirements.txt
python -m pytest -q
python experiments/budget/get_openmrg.py        # 下载到data/OpenMRG并校验MD5
python experiments/budget/prepare_openmrg.py    # 生成分钟缓存
python experiments/budget/analyze_openmrg.py    # 生成分析表和offset.json
python experiments/budget/pilot_reference.py    # 运行NR DM-RS时域仿真
python experiments/budget/summarize.py
python experiments/budget/make_figures.py       # 生成paper/figures_zh和图例记录
python experiments/budget/summarize.py          # 将图1事件信息纳入summary.json
python paper/build_manuscript.py --lang zh      # 生成中文Word和PDF
```

随机种子、筛选阈值、参考窗口、雨强分档、月份划分和路径分组设置集中在`experiments/budget/config.json`。重新生成结果后，应再次执行归档核验命令，并同步更新文稿中的来源与重算状态。

## 中文文稿构建

文稿源文件为`paper/content_zh.md`，数值占位符从`results/budget/summary.json`读取。中文稿沿用`paper/template/reference.docx`的页设置，运行字体按下表配置。

| 用途 | 字体 |
|---|---|
| 中文正文 | AR PL UMing CN |
| 中文标题 | Noto Sans CJK SC |
| 拉丁文字 | Nimbus Roman |
| 数学公式 | Latin Modern Math |

构建环境需要XeLaTeX及`unicode-math`、`amsmath`、`standalone`宏包，Poppler的`pdfinfo`和`pdftoppm`，以及提供`soffice`的LibreOffice。修订后的构建器兼容稿件所用的行内数学宏，每次PDF转换使用独立的LibreOffice用户配置目录。

```bash
python paper/build_manuscript.py --lang zh
```

生成文件为`paper/Manuscript_zh.docx`和`paper/Manuscript_zh.pdf`。复用归档图表构建的PDF会保留量化历史结果待重算的标记，完成完整流程后再更新相应状态。

## 数据来源

OpenMRG由Andersson等发表于2022年《Earth System Science Data》，以CC BY-SA 4.0许可发布。下载入口与数据版本见[OpenMRG数据记录](https://doi.org/10.5281/zenodo.7107689)。原始数据存放于`data/OpenMRG`，仓库中的`results/budget/`提供派生结果和核验记录。

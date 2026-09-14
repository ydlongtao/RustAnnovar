# 兼容性状态

本文件区分已经实现并经过本地注册版 ANNOVAR 测试的行为，以及开放测试阶段仍需补齐的兼容细节。

> **开放测试状态：** 当前版本用于研究验证和性能测试。复杂或高风险结果应使用原版 ANNOVAR 或其他成熟工具复核。

## 已实现

- VCF 4.x 文本及 gzip 输入；多等位记录拆分；AVinput 读写。
- 内部零起点半开区间；AVinput 边界转换；常见 VCF Indel 锚碱基移除。
- ANNOVAR generic filter 数据库的 Chr/Start/End/Ref/Alt 精确匹配。
- generic、BED 风格及 UCSC bin 风格区间数据库的重叠查询。
- UCSC refGene 风格基因模型的外显子、内含子、UTR、剪接、上下游分类。
- 有转录本 FASTA 时的编码 SNV 同义、非同义、stopgain、stoploss 计算。
- 多协议 TSV/CSV 汇总和 VCF INFO 注释；输入附加列保留。
- 数据库来源指纹、陈旧检测、序列提取和表格等值筛选。

## 原版基线

本地基准固定为 `annotate_variation.pl` 2025-03-02、`table_annovar.pl` 2022-08-02。脚本 SHA-256 与官方 `ex1` refGene 输出保存在 `tests/fixtures/perl_expected/`。忽略测试会加载注册安装包，核对坐标、功能分类、基因、全部外显子 SNV 的转录本和蛋白变化，以及原包 generic 过滤与 GFF3 区域示例。

## 尚待完整兼容

- ANNOVAR 对复杂多外显子 Indel、HGVS 边界和异常 ORF 的全部历史细节；当前真实基线已确认 SNV 核心后果一致。
- 各个特殊 `dbtype` 的隐式列规则、阈值规则及遗留格式。
- `convert2annovar.pl` 支持的 VCF 之外的历史测序格式。
- `coding_change.pl` 的完整突变蛋白 FASTA 输出格式。
- ANNOVAR 下载服务的注册、授权及完整数据库目录协议。
- 命令行参数逐项兼容；本项目当前保证的是主要数据与结果字段兼容。

## 性能边界

注释行在多个 CPU 核心并行执行，结果顺序稳定。`db index --kind filter` 会按 1 Mb 基因组区块记录源文件字节范围，查询时只读取输入变异涉及的区块；没有索引、索引陈旧或 gzip 数据库会安全回退到完整加载。200 GB 级数据库仍需用真实数据测量索引体积、冷缓存吞吐量和并发 I/O，才满足正式 WGS 发布条件。

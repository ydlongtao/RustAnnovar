# RustAnnovar

[English](README.md) | **简体中文**

[![CI](https://github.com/ydlongtao/RustAnnovar/actions/workflows/ci.yml/badge.svg)](https://github.com/ydlongtao/RustAnnovar/actions/workflows/ci.yml)
[![Gem Version](https://badge.fury.io/rb/rust-annovar.svg)](https://rubygems.org/gems/rust-annovar)
[![Open Beta](https://img.shields.io/badge/status-open%20beta-orange)](https://github.com/ydlongtao/RustAnnovar/issues)
[![Rust](https://img.shields.io/badge/Rust-1.85%2B-000000?logo=rust)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

**A high-performance ANNOVAR-compatible variant annotation engine written in Rust.**

`RustAnnovar` 是一个使用 Rust 编写的高性能基因组变异注释工具。它直接读取现有 ANNOVAR `humandb` 数据文件，支持 VCF/AVinput 转换、基因注释、区域注释、过滤注释、多数据库汇总和 VCF 回写，并提供面向大型数据库的索引能力。

> [!WARNING]
> **本软件正在开放测试（Open Beta）。** 当前版本适合功能验证、性能测试和非关键研究流程。SNV 核心后果、generic filter 和 GFF3 区域注释已经与本地 ANNOVAR 基线核对；复杂 Indel、完整 HGVS、ncRNA 分类细节和部分历史数据库协议仍在完善。请在研究或临床决策前用原版 ANNOVAR 或其他成熟工具复核结果，并通过 [Issues](https://github.com/ydlongtao/RustAnnovar/issues) 报告差异。

## 为什么使用 RustAnnovar

- **原生 Rust 引擎**：核心注释流程不调用 Perl，提供内存安全和稳定并行执行。
- **复用 humandb**：按 `hg19_refGene.txt`、`hg38_clinvar.txt` 等 ANNOVAR 命名方式发现数据库。
- **三类核心查询**：精确变异匹配、基因组区间重叠、转录本与编码后果计算。
- **常用输入输出**：读取 VCF、gzip VCF 和 AVinput，输出 TSV、CSV 或带 INFO 注释的 VCF。
- **大型数据库索引**：过滤数据库可建立 1 Mb 分块索引，只读取输入位点涉及的数据区块。
- **可审计兼容性**：仓库包含自动化测试和与注册版 ANNOVAR 对照的回归测试框架。

## 运行速度

在 Apple M1 8 核、hg19 refGene、热文件缓存、release 构建条件下，多次运行取墙钟时间中位数：

| 场景 | 原版 ANNOVAR（Perl） | RustAnnovar（Rust） | 加速比 |
|---|---:|---:|---:|
| 13 条已核对 SNV | 2.40 s | 0.29 s | **8.28×** |
| 26,000 行 SNV，默认设置 | 4.68 s | 0.37 s | **12.65×** |
| 26,000 行 SNV，单线程 | 4.70 s | 0.47 s | **10.00×** |
| 21 个变异查询 25,688 个 GFF3 区域 | 0.14 s | 0.01 s | **约 14×** |

26,000 行测试中，最大常驻内存由 400 MiB 降至 360.3 MiB，减少约 9.9%。该输入由 13 条已经核对的 SNV 重复构成，用于稳定测量启动、注释和输出成本；它不代表真实 WES 的独立位点分布。完整环境、原始计时和兼容性边界见 [性能报告](docs/BENCHMARK_2026-09-14.md)。

## 安装

### 环境要求

- Linux 或 macOS；Windows 可从源码构建
- Rust 1.85 或更高版本
- 用户自行取得的 ANNOVAR `humandb` 数据库；演示命令不需要注册数据库

没有 Rust 时，先安装官方工具链：

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

### 从 GitHub 安装

```bash
cargo install --git https://github.com/ydlongtao/RustAnnovar.git --locked
rust-annovar --version
```

### 从 RubyGems 安装

[RubyGems 项目页](https://rubygems.org/gems/rust-annovar) · [0.1.0.beta.1 版本页](https://rubygems.org/gems/rust-annovar/versions/0.1.0.beta.1) · [直接下载 .gem 安装包](https://rubygems.org/gems/rust-annovar-0.1.0.beta.1.gem)

RubyGems 包会在安装时使用 Cargo 编译 Rust 可执行文件，因此仍需预先安装 Rust 工具链和本机链接器。Ruby 只负责封装和启动编译后的程序。

```bash
gem install rust-annovar --pre
rust-annovar --version
```

RubyGems 版本号显示为 `0.1.0.beta.1`，命令行版本号显示为 `0.1.0-beta.1`。

安装指定测试版本：

```bash
gem install rust-annovar --version 0.1.0.beta.1 --pre
```

也可以使用上方链接下载安装包，再执行：

```bash
gem install ./rust-annovar-0.1.0.beta.1.gem
rust-annovar --help
```

下载的安装包仍需 Cargo 和本机链接器进行编译，安装时 Cargo 可能联网获取 Rust 依赖。包内不包含 ANNOVAR 数据库；可以使用下方合成示例体验，或自行提供兼容数据库。

### 克隆源码构建

```bash
git clone https://github.com/ydlongtao/RustAnnovar.git
cd RustAnnovar
cargo build --release --locked
cargo test --all-features
./target/release/rust-annovar --version
```

## 五分钟快速体验

仓库包含一个完全合成的小型 VCF 和过滤数据库：

```bash
rust-annovar table \
  examples/demo.vcf examples/humandb \
  --build hg38 \
  --protocol demo \
  --operation f \
  --vcf-input \
  --output demo.multianno.tsv \
  --vcf-output demo.annotated.vcf

cat demo.multianno.tsv
```

预期第一条变异命中 `Pathogenic`，第二条变异显示缺失值 `.`。

## 使用现有 humandb

目录按 ANNOVAR 的文件名规则组织：

```text
humandb/
├── hg38_refGene.txt
├── hg38_refGeneMrna.fa
├── hg38_cytoBand.txt
└── hg38_clinvar.txt
```

组合基因、区域和过滤注释：

```bash
rust-annovar table sample.vcf humandb \
  --build hg38 \
  --protocol refGene,cytoBand,clinvar \
  --operation g,r,f \
  --vcf-input \
  --output sample.hg38_multianno.tsv \
  --vcf-output sample.hg38_multianno.vcf
```

`--protocol` 与 `--operation` 必须一一对应；`g`、`r`、`f` 分别代表 gene、region 和 filter。输出列按照协议请求顺序排列，VCF 输出保留原始样本与 FORMAT 字段。

## 常用命令

### VCF 转换为 AVinput

```bash
rust-annovar convert sample.vcf.gz \
  --include-info \
  --output sample.avinput
```

程序会拆分多等位记录，并处理常见 VCF Indel 的锚碱基。

### 单个过滤数据库

```bash
rust-annovar annotate sample.avinput humandb/hg38_clinvar.txt \
  --operation filter \
  --protocol clinvar \
  --output sample.clinvar.tsv
```

过滤注释要求染色体、坐标、Ref 和 Alt 精确匹配。

### 单个区域数据库

```bash
rust-annovar annotate sample.avinput humandb/hg38_cytoBand.txt \
  --operation region \
  --protocol cytoBand \
  --output sample.cytoband.tsv
```

### refGene 基因注释

```bash
rust-annovar annotate sample.avinput humandb/hg38_refGene.txt \
  --operation gene \
  --protocol refGene \
  --fasta humandb/hg38_refGeneMrna.fa \
  --output sample.refgene.tsv
```

### 为大型过滤数据库建立索引

```bash
rust-annovar db index humandb/hg38_dbnsfp.txt --kind filter
rust-annovar db check humandb/hg38_dbnsfp.fai.json
```

源数据库变更后，`db check` 会报告索引过期。gzip 数据库会安全回退到完整加载。

### 提取序列和筛选结果

```bash
rust-annovar sequence regions.avinput reference.fa --output regions.fa

rust-annovar reduce sample.hg38_multianno.tsv \
  --column Func.refGene \
  --equals exonic \
  --output sample.exonic.tsv
```

运行 `rust-annovar <子命令> --help` 可查看完整参数。

## 子命令概览

| 子命令 | 功能 |
|---|---|
| `convert` | VCF/gzip VCF 转换为 AVinput |
| `annotate` | 执行单个 gene、region 或 filter 注释 |
| `table` | 组合多个数据库并输出 TSV、CSV 或 VCF |
| `db` | 建立和校验索引、列出或下载公开数据库 |
| `sequence` | 按 AVinput 区间提取 FASTA 序列 |
| `coding-change` | 输出基因与编码后果注释 |
| `reduce` | 按指定注释列筛选多注释表格 |

## 测试与兼容性

```bash
cargo fmt -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
```

如果本机已有注册版 ANNOVAR，可运行对照测试：

```bash
ANNOVAR_HOME=/path/to/annovar \
  cargo test --test annovar_compat -- --ignored
```

已实现范围和已知差异见 [兼容性状态](docs/COMPATIBILITY.md)。发现差异时，请在 issue 中同时提供最小化输入、数据库版本、构建版本、原版命令和两份输出；请勿上传受许可限制的数据库或可识别个体的基因组数据。

## ANNOVAR 与数据库说明

本项目是独立实现，不隶属于 ANNOVAR，也不分发 ANNOVAR Perl 程序或注册数据库。ANNOVAR 及相关数据库可能有各自的学术或商业许可；用户需要自行确认使用资格。项目名称用于说明数据格式和结果兼容目标。

## 许可证

源代码采用 [MIT](LICENSE-MIT) 或 [Apache-2.0](LICENSE-APACHE) 双许可证。用户可任选其一。

## 致谢

项目结构和公开文档风格参考了 [Huang-lab/fastVEP](https://github.com/Huang-lab/fastVEP)。感谢 ANNOVAR 作者和变异注释社区建立的数据格式、数据库与验证案例。

# ANNOVAR 对照基线

此目录用于保存固定版本 ANNOVAR 产生的期望结果。官方安装包需要注册获取，因此仓库不分发其 Perl 源码。

取得安装包后，设置 `ANNOVAR_HOME`，再运行：

```bash
scripts/capture_perl_baseline.sh
```

脚本会记录版本清单和 SHA-256，并把官方示例输出写入本目录。结果进入版本控制前，应移除日期、绝对路径等不稳定字段。

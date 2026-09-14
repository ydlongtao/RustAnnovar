# 参与开放测试

感谢测试 RustAnnovar。当前最需要的是能够复现的 ANNOVAR 兼容性报告，以及真实 WES/WGS 工作负载下的性能数据。

提交 issue 时请包含：

1. `rust-annovar --version`、操作系统和 CPU。
2. hg19、hg38 或其他构建版本，以及数据库文件名和版本。
3. 最小化、去标识化的 VCF 或 AVinput。
4. RustAnnovar 与原版 ANNOVAR 的完整命令。
5. 期望输出、实际输出和具体差异字段。

请勿提交 ANNOVAR Perl 程序、受许可限制的数据库、访问令牌或可识别个体的基因组数据。

代码变更应通过：

```bash
cargo fmt -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
```

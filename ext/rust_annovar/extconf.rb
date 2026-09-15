# frozen_string_literal: true

require "rbconfig"
require "shellwords"

manifest = File.expand_path("../../Cargo.toml", __dir__)
cargo = ENV.fetch("CARGO", "cargo")
unless system(cargo, "--version", out: File::NULL)
  abort "Cargo is required to install rust-annovar. Install Rust from https://rustup.rs/"
end

manifest_arg = Shellwords.escape(manifest)
cargo_arg = Shellwords.escape(cargo)
executable = "rust-annovar#{RbConfig::CONFIG.fetch("EXEEXT", "")}"
target_dir = File.expand_path("target", __dir__)
install_dir = File.join(RbConfig::CONFIG.fetch("sitearchdir"), "rust_annovar")
ruby_arg = Shellwords.escape(RbConfig.ruby)
target_arg = Shellwords.escape(target_dir)
install_arg = Shellwords.escape(install_dir)
binary_arg = Shellwords.escape(File.join(target_dir, "release", executable))

File.write(
  "Makefile",
  <<~MAKEFILE
    all:
    \tCARGO_TARGET_DIR=#{target_arg} #{cargo_arg} build --release --locked --bin rust-annovar --manifest-path #{manifest_arg}

    install:
    \t#{ruby_arg} -rfileutils -e 'FileUtils.mkdir_p(ARGV.fetch(0)); FileUtils.install(ARGV.fetch(1), ARGV.fetch(0), mode: 0755)' #{install_arg} #{binary_arg}

    clean:
    \t#{ruby_arg} -rfileutils -e 'FileUtils.rm_rf(ARGV.fetch(0))' #{target_arg}
  MAKEFILE
)

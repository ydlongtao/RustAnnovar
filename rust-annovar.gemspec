# frozen_string_literal: true

require_relative "lib/rust_annovar/version"

Gem::Specification.new do |spec|
  spec.name = "rust-annovar"
  spec.version = RustAnnovar::VERSION
  spec.authors = ["Longtao Huangfu"]
  spec.email = ["ydlongtao@gmail.com"]

  spec.summary = "Rust implementation of core ANNOVAR-compatible workflows"
  spec.description = <<~DESCRIPTION.strip
    RustAnnovar is an open-beta genomic variant annotation command-line tool.
    This source gem compiles the Rust executable during installation and reads
    supported ANNOVAR humandb formats supplied separately by the user.
  DESCRIPTION
  spec.homepage = "https://github.com/ydlongtao/RustAnnovar"
  spec.licenses = ["MIT", "Apache-2.0"]
  spec.required_ruby_version = Gem::Requirement.new(">= 2.6")

  spec.metadata = {
    "github_repo" => "ssh://github.com/ydlongtao/RustAnnovar",
    "source_code_uri" => "https://github.com/ydlongtao/RustAnnovar",
    "documentation_uri" => "https://github.com/ydlongtao/RustAnnovar#readme",
    "bug_tracker_uri" => "https://github.com/ydlongtao/RustAnnovar/issues"
  }

  spec.files = Dir[
    "bin/rust-annovar",
    "ext/rust_annovar/extconf.rb",
    "lib/**/*.rb",
    "src/**/*.rs",
    "crates/**/Cargo.toml",
    "crates/**/src/**/*.rs",
    "Cargo.toml",
    "Cargo.lock",
    "README.md",
    "ARCHITECTURE.md",
    "PERFORMANCE.md",
    "CHANGELOG.md",
    "docs/CADD*.md",
    "scripts/*cadd*",
    "LICENSE-MIT",
    "LICENSE-APACHE"
  ]
  spec.bindir = "bin"
  spec.executables = ["rust-annovar"]
  spec.require_paths = ["lib"]
  spec.extensions = ["ext/rust_annovar/extconf.rb"]
end

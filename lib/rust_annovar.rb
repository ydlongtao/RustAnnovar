# frozen_string_literal: true

require "rbconfig"
require_relative "rust_annovar/version"

module RustAnnovar
  class Error < StandardError; end

  def self.executable
    specification = Gem.loaded_specs.fetch("rust-annovar")
    suffix = RbConfig::CONFIG.fetch("EXEEXT", "")
    executable_name = "rust-annovar#{suffix}"
    candidates = [
      File.join(specification.extension_dir, "rust_annovar", executable_name),
      File.join(RbConfig::CONFIG.fetch("sitearchdir"), "rust_annovar", executable_name),
      File.join(specification.full_gem_path, "ext", "rust_annovar", "target", "release", executable_name)
    ]
    candidates.concat(
      Dir.glob(
        File.join(
          specification.base_dir,
          "extensions",
          "**",
          specification.full_name,
          "rust_annovar",
          executable_name
        )
      )
    )
    path = candidates.find { |candidate| File.file?(candidate) && File.executable?(candidate) }
    return path if path

    raise Error, "compiled rust-annovar executable was not found in the RubyGems extension directories"
  end
end

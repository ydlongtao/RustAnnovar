# frozen_string_literal: true

require "rbconfig"
require_relative "rust_annovar/version"

module RustAnnovar
  class Error < StandardError; end

  def self.executable
    specification = Gem.loaded_specs.fetch("rust-annovar")
    suffix = RbConfig::CONFIG.fetch("EXEEXT", "")
    path = File.join(specification.extension_dir, "rust_annovar", "rust-annovar#{suffix}")
    return path if File.file?(path) && File.executable?(path)

    raise Error, "compiled rust-annovar executable was not found at #{path}"
  end
end

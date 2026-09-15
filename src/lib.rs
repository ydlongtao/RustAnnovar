//! Core annotation engine for `RustAnnovar`.
//!
//! Coordinates are always represented internally as zero-based half-open
//! intervals. Parsers and writers are responsible for conversion at the edge.

pub mod database;
pub mod gene;
pub mod io;
pub mod model;
pub mod pipeline;

pub use model::{Annotation, AnnotationKind, Variant};

pub mod batch;
pub mod disk_index;
pub mod interval;
pub mod reference;
pub mod stream;

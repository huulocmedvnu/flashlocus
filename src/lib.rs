//! flashLocus: summary-statistic fine-mapping (exact port of susieR 0.14.2 susie_rss), colocalisation (exact
//! port of coloc 5.2.3 coloc.bf_bf / coloc.susie) and exact LD from PLINK files, in Rust.
pub mod atlas;
pub mod bfile;
pub mod brent;
pub mod coloc;
pub mod compare;
pub mod ldmat;
pub mod prep;
pub mod rng;
pub mod sim;
pub mod simd;
pub mod susie;
pub mod topo;
pub mod cli;
#[cfg(feature = "python")]
mod python;

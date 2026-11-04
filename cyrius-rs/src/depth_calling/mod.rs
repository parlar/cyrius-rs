pub mod bin_count;
pub mod copy_number_call;
pub mod gmm;
pub mod haplotype;
pub mod snp_count;
pub mod utilities;

/// Minimum base quality for an allele-discriminating base: the
/// `bcftools mpileup` default (`--min-BQ 13`). A pattern counts only if every
/// base it was matched on reaches it. Used by both the fragment-phasing path
/// (`caller::phase_disambiguate`) and the targeted-rule path
/// (`depth_calling::haplotype`), so it lives here, in the lower layer.
pub const PHASE_MIN_BASE_QUALITY: u8 = 13;

/// The value SAM/BAM stores when base qualities are **omitted**. It is not a
/// quality of 255: it means there is no quality information, so it can never
/// satisfy `PHASE_MIN_BASE_QUALITY` (review #5).
pub const QUAL_UNAVAILABLE: u8 = 0xFF;

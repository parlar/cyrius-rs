/// Embedded data files for BCyrius.
/// All data is compiled into the binary — no external files needed at runtime.

pub const STAR_TABLE: &str = include_str!("../data/star_table.txt");
pub const GMM_PARAMS: &str = include_str!("../data/CYP2D6_gmm.txt");
pub const REGION_BED: &str = include_str!("../data/CYP2D6_region_38.bed");
pub const SNP_FILE: &str = include_str!("../data/CYP2D6_SNP_38.txt");
pub const TARGET_VARIANT: &str = include_str!("../data/CYP2D6_target_variant_38.txt");
pub const TARGET_VARIANT_HOMO: &str =
    include_str!("../data/CYP2D6_target_variant_homology_region_38.txt");
pub const HAPLOTYPE_FILE: &str = include_str!("../data/CYP2D6_haplotype_38.txt");
pub const HAPLOTYPE_FUNC: &str = include_str!("../data/CYP2D6_haplotypes_functionality.txt");
pub const FREQ_HAP: &str = include_str!("../data/CYP2D6_frequency_table_haplotypes.txt");
pub const FREQ_DIP: &str = include_str!("../data/CYP2D6_frequency_table_diplotypes.txt");

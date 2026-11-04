/// The `Filter` value for a final call, from the first applicable reason it is not
/// a clean single answer. A masked set that a later step narrowed to one answer
/// (`unknown_sites_ambiguous`) reads `Unknown_sites_resolved`, never `PASS`.
pub fn genotype_filter(
    call: &str,
    high_cn_low_confidence: bool,
    unknown_sites_unmatched: bool,
    unknown_sites_ambiguous: bool,
    is_fuzzy: bool,
) -> &'static str {
    if call.contains(';') {
        "More_than_one_possible_genotype"
    } else if !call.contains('/') {
        "Not_assigned_to_haplotypes"
    } else if high_cn_low_confidence {
        "LowQ_high_CN"
    } else if unknown_sites_unmatched {
        "Unknown_sites_unmatched"
    } else if unknown_sites_ambiguous {
        "Unknown_sites_resolved"
    } else if is_fuzzy {
        "Fuzzy_match"
    } else {
        "PASS"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_branch_in_order() {
        assert_eq!(genotype_filter("*1/*2;*1/*3", true, true, true, true), "More_than_one_possible_genotype");
        assert_eq!(genotype_filter("*1", true, true, true, true), "Not_assigned_to_haplotypes");
        assert_eq!(genotype_filter("*1/*2", true, true, true, true), "LowQ_high_CN");
        assert_eq!(genotype_filter("*1/*2", false, true, true, true), "Unknown_sites_unmatched");
        assert_eq!(genotype_filter("*1/*2", false, false, true, true), "Unknown_sites_resolved");
        assert_eq!(genotype_filter("*1/*2", false, false, false, true), "Fuzzy_match");
        assert_eq!(genotype_filter("*1/*2", false, false, false, false), "PASS");
    }

    #[test]
    fn either_unknown_sites_flag_is_never_pass() {
        for call in ["*1/*2", "*1/*2;*1/*3", "*1"] {
            for bits in 0..16u8 {
                let (hq, un, am, fz) = (bits & 1 != 0, bits & 2 != 0, bits & 4 != 0, bits & 8 != 0);
                if un || am {
                    assert_ne!(genotype_filter(call, hq, un, am, fz), "PASS", "{call} {bits:04b}");
                }
            }
        }
    }
}

use crate::types::StarCombinations;
use std::collections::HashMap;

/// Exon 9 gene conversion alleles
const EXON9GC_ALLELES: &[&str] = &["*36", "*4.013", "*57", "*83", "*141"];

fn exon9gc_pair_alleles() -> HashMap<&'static str, &'static str> {
    let mut m = HashMap::new();
    m.insert("*36", "*10");
    m.insert("*4.013", "*4");
    m
}

/// Update variant-to-diplotype dictionary.
fn make_hap_dic(variant_list: &[&str], star_set: &str, hap_dic: &mut HashMap<String, Vec<String>>) {
    let filtered: Vec<&str> = variant_list.iter().copied().filter(|&v| v != "NA").collect();
    let var_list_joined = if filtered.is_empty() {
        "NA".to_string()
    } else {
        filtered.join("_")
    };
    let entry = hap_dic.entry(var_list_joined).or_insert_with(Vec::new);
    if !entry.contains(&star_set.to_string()) {
        entry.push(star_set.to_string());
    }
}

/// Return all the possible variant tables based on the star allele definition file.
pub fn get_hap_table(star_table_content: &str) -> StarCombinations {
    // CN=1, one copy
    let mut dhap: HashMap<String, String> = HashMap::new();
    let mut dstar: HashMap<String, String> = HashMap::new();

    for line in star_table_content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let star_id = fields[0];
        let mut variant_list: Vec<&str> = fields[1..fields.len() - 1].to_vec();
        variant_list.sort();
        let var_list_joined = variant_list.join("_");
        dhap.entry(var_list_joined.clone())
            .or_insert_with(|| star_id.to_string());
        dstar
            .entry(star_id.to_string())
            .or_insert_with(|| var_list_joined);
    }

    let star_keys: Vec<String> = dstar.keys().cloned().collect();

    // CN=2, two copies
    let mut dhap2: HashMap<String, Vec<String>> = HashMap::new();
    for star1 in &star_keys {
        for star2 in &star_keys {
            let mut variant_list: Vec<&str> = dstar[star1]
                .split('_')
                .chain(dstar[star2].split('_'))
                .collect();
            variant_list.sort();
            let star_set = {
                let mut s = vec![star1.as_str(), star2.as_str()];
                s.sort();
                s.join("_")
            };
            make_hap_dic(&variant_list, &star_set, &mut dhap2);
        }
    }

    // CN=3, three copies (for exon9hyb cases)
    let mut dhap3: HashMap<String, Vec<String>> = HashMap::new();
    for star1 in &star_keys {
        for star2 in &star_keys {
            for star3 in &star_keys {
                let mut variant_list: Vec<&str> = dstar[star1]
                    .split('_')
                    .chain(dstar[star2].split('_'))
                    .chain(dstar[star3].split('_'))
                    .collect();
                variant_list.sort();
                let star_set = {
                    let mut s = vec![star1.as_str(), star2.as_str(), star3.as_str()];
                    s.sort();
                    s.join("_")
                };
                make_hap_dic(&variant_list, &star_set, &mut dhap3);
            }
        }
    }

    // CN=3, three copies, no hybrid, only duplication (duplicated copies identical)
    let mut dhap3pair: HashMap<String, Vec<String>> = HashMap::new();
    for star1 in &star_keys {
        for star2 in &star_keys {
            let mut variant_list: Vec<&str> = dstar[star1]
                .split('_')
                .chain(dstar[star2].split('_'))
                .chain(dstar[star2].split('_'))
                .collect();
            variant_list.sort();
            let star_set = {
                let mut s = vec![star1.as_str(), star2.as_str(), star2.as_str()];
                s.sort();
                s.join("_")
            };
            make_hap_dic(&variant_list, &star_set, &mut dhap3pair);
        }
    }

    // CN=4, four copies (star1x2 + star2x2 or star1 + star2x3)
    let mut dhap4pair: HashMap<String, Vec<String>> = HashMap::new();
    for star1 in &star_keys {
        for star2 in &star_keys {
            // star1x2 + star2x2
            let mut vl: Vec<&str> = dstar[star1]
                .split('_')
                .chain(dstar[star1].split('_'))
                .chain(dstar[star2].split('_'))
                .chain(dstar[star2].split('_'))
                .collect();
            vl.sort();
            let ss = {
                let mut s = vec![star1.as_str(), star1.as_str(), star2.as_str(), star2.as_str()];
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap4pair);

            // star1 + star2x3
            let mut vl: Vec<&str> = dstar[star1]
                .split('_')
                .chain(dstar[star2].split('_'))
                .chain(dstar[star2].split('_'))
                .chain(dstar[star2].split('_'))
                .collect();
            vl.sort();
            let ss = {
                let mut s = vec![star1.as_str(), star2.as_str(), star2.as_str(), star2.as_str()];
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap4pair);
        }
    }

    // CN=5
    let mut dhap5pair: HashMap<String, Vec<String>> = HashMap::new();
    for star1 in &star_keys {
        for star2 in &star_keys {
            // star1x2 + star2x3
            let mut vl: Vec<&str> = std::iter::repeat(dstar[star1].as_str())
                .take(2)
                .flat_map(|s| s.split('_'))
                .chain(std::iter::repeat(dstar[star2].as_str()).take(3).flat_map(|s| s.split('_')))
                .collect();
            vl.sort();
            let ss = {
                let mut s: Vec<&str> = std::iter::repeat(star1.as_str())
                    .take(2)
                    .chain(std::iter::repeat(star2.as_str()).take(3))
                    .collect();
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap5pair);

            // star1 + star2x4
            let mut vl: Vec<&str> = dstar[star1]
                .split('_')
                .chain(std::iter::repeat(dstar[star2].as_str()).take(4).flat_map(|s| s.split('_')))
                .collect();
            vl.sort();
            let ss = {
                let mut s: Vec<&str> = std::iter::once(star1.as_str())
                    .chain(std::iter::repeat(star2.as_str()).take(4))
                    .collect();
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap5pair);
        }
    }

    // CN=6
    let mut dhap6pair: HashMap<String, Vec<String>> = HashMap::new();
    for star1 in &star_keys {
        for star2 in &star_keys {
            // star1x2 + star2x4
            let mut vl: Vec<&str> = std::iter::repeat(dstar[star1].as_str())
                .take(2)
                .flat_map(|s| s.split('_'))
                .chain(std::iter::repeat(dstar[star2].as_str()).take(4).flat_map(|s| s.split('_')))
                .collect();
            vl.sort();
            let ss = {
                let mut s: Vec<&str> = std::iter::repeat(star1.as_str())
                    .take(2)
                    .chain(std::iter::repeat(star2.as_str()).take(4))
                    .collect();
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap6pair);

            // star1 + star2x5
            let mut vl: Vec<&str> = dstar[star1]
                .split('_')
                .chain(std::iter::repeat(dstar[star2].as_str()).take(5).flat_map(|s| s.split('_')))
                .collect();
            vl.sort();
            let ss = {
                let mut s: Vec<&str> = std::iter::once(star1.as_str())
                    .chain(std::iter::repeat(star2.as_str()).take(5))
                    .collect();
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap6pair);

            // star1x3 + star2x3
            let mut vl: Vec<&str> = std::iter::repeat(dstar[star1].as_str())
                .take(3)
                .flat_map(|s| s.split('_'))
                .chain(std::iter::repeat(dstar[star2].as_str()).take(3).flat_map(|s| s.split('_')))
                .collect();
            vl.sort();
            let ss = {
                let mut s: Vec<&str> = std::iter::repeat(star1.as_str())
                    .take(3)
                    .chain(std::iter::repeat(star2.as_str()).take(3))
                    .collect();
                s.sort();
                s.join("_")
            };
            make_hap_dic(&vl, &ss, &mut dhap6pair);
        }
    }

    // exon9hybx2: limit search to EXON9GC_ALLELES
    let mut dhap_exon9_x2: HashMap<String, Vec<String>> = HashMap::new();
    for &e1 in EXON9GC_ALLELES {
        for &e2 in EXON9GC_ALLELES {
            if dstar.contains_key(e1) && dstar.contains_key(e2) {
                for star3 in &star_keys {
                    for star4 in &star_keys {
                        let mut vl: Vec<&str> = dstar[e1]
                            .split('_')
                            .chain(dstar[e2].split('_'))
                            .chain(dstar[star3].split('_'))
                            .chain(dstar[star4].split('_'))
                            .collect();
                        vl.sort();
                        let ss = {
                            let mut s = vec![e1, e2, star3.as_str(), star4.as_str()];
                            s.sort();
                            s.join("_")
                        };
                        make_hap_dic(&vl, &ss, &mut dhap_exon9_x2);
                    }
                }
            }
        }
    }

    // exon9hybx3: limit search to EXON9GC_ALLELES
    let mut dhap_exon9_x3: HashMap<String, Vec<String>> = HashMap::new();
    for &e1 in EXON9GC_ALLELES {
        for &e2 in EXON9GC_ALLELES {
            for &e3 in EXON9GC_ALLELES {
                if dstar.contains_key(e1) && dstar.contains_key(e2) && dstar.contains_key(e3) {
                    for star4 in &star_keys {
                        for star5 in &star_keys {
                            let mut vl: Vec<&str> = dstar[e1]
                                .split('_')
                                .chain(dstar[e2].split('_'))
                                .chain(dstar[e3].split('_'))
                                .chain(dstar[star4].split('_'))
                                .chain(dstar[star5].split('_'))
                                .collect();
                            vl.sort();
                            let ss = {
                                let mut s =
                                    vec![e1, e2, e3, star4.as_str(), star5.as_str()];
                                s.sort();
                                s.join("_")
                            };
                            make_hap_dic(&vl, &ss, &mut dhap_exon9_x3);
                        }
                    }
                }
            }
        }
    }

    // exon9hybx4: limit to *10 and *36
    let mut dhap_exon9_x4: HashMap<String, Vec<String>> = HashMap::new();
    if dstar.contains_key("*10") && dstar.contains_key("*36") {
        let mut vl: Vec<&str> = std::iter::repeat(dstar["*10"].as_str())
            .take(2)
            .flat_map(|s| s.split('_'))
            .chain(
                std::iter::repeat(dstar["*36"].as_str())
                    .take(4)
                    .flat_map(|s| s.split('_')),
            )
            .collect();
        vl.sort();
        let ss = {
            let mut s = vec!["*10", "*10", "*36", "*36", "*36", "*36"];
            s.sort();
            s.join("_")
        };
        make_hap_dic(&vl, &ss, &mut dhap_exon9_x4);
    }

    // dup_exon9hyb: for exon9hyb part, limit to EXON9GC_PAIR_ALLELES
    let mut dhap_dup_exon9: HashMap<String, Vec<String>> = HashMap::new();
    let pair_alleles = exon9gc_pair_alleles();
    for (&e1, &e2) in &pair_alleles {
        if dstar.contains_key(e1) && dstar.contains_key(e2) {
            for star3 in &star_keys {
                for star4 in &star_keys {
                    let mut vl: Vec<&str> = dstar[e1]
                        .split('_')
                        .chain(dstar[e2].split('_'))
                        .chain(dstar[star3].split('_'))
                        .chain(dstar[star4].split('_'))
                        .collect();
                    vl.sort();
                    let ss = {
                        let mut s = vec![e1, e2, star3.as_str(), star4.as_str()];
                        s.sort();
                        s.join("_")
                    };
                    make_hap_dic(&vl, &ss, &mut dhap_dup_exon9);
                }
            }
        }
    }

    StarCombinations {
        dhap,
        dhap2,
        dhap3,
        dhap3pair,
        dhap4pair,
        dhap5pair,
        dhap6pair,
        dhap_exon9_x2,
        dhap_exon9_x3,
        dhap_exon9_x4,
        dhap_dup_exon9,
    }
}

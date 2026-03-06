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

/// Merge variant ID slices, sort by ID (which equals alphabetical order),
/// convert to "_"-joined string key. IDs must be assigned in alphabetical order.
fn merge_to_key(slices: &[&[u16]], id_to_var: &[&str]) -> String {
    let total: usize = slices.iter().map(|s| s.len()).sum();
    let mut merged = Vec::with_capacity(total);
    for s in slices {
        merged.extend_from_slice(s);
    }
    merged.sort_unstable();

    // Filter NA (which has a known ID) and build key
    let mut first = true;
    let mut key = String::with_capacity(total * 20);
    for &id in &merged {
        let name = id_to_var[id as usize];
        if name == "NA" {
            continue;
        }
        if !first {
            key.push('_');
        }
        key.push_str(name);
        first = false;
    }
    if key.is_empty() {
        "NA".to_string()
    } else {
        key
    }
}

/// Build a sorted star set string from star names.
fn make_star_set(stars: &[&str]) -> String {
    let mut s: Vec<&str> = stars.to_vec();
    s.sort();
    s.join("_")
}

fn insert_hap(key: &str, star_set: &str, hap_dic: &mut HashMap<String, Vec<String>>) {
    let entry = hap_dic.entry(key.to_string()).or_insert_with(Vec::new);
    if !entry.contains(&star_set.to_string()) {
        entry.push(star_set.to_string());
    }
}

/// Return all the possible variant tables based on the star allele definition file.
pub fn get_hap_table(star_table_content: &str) -> StarCombinations {
    // Phase 1: Parse and collect all unique variant names
    let mut all_var_names: Vec<String> = Vec::new();
    let mut var_name_set: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut star_entries: Vec<(&str, Vec<&str>)> = Vec::new();

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
        let variant_list: Vec<&str> = fields[1..fields.len() - 1].to_vec();
        for &v in &variant_list {
            if var_name_set.insert(v.to_string()) {
                all_var_names.push(v.to_string());
            }
        }
        star_entries.push((star_id, variant_list));
    }

    // Sort variant names alphabetically and assign IDs in that order
    // This ensures: id_a < id_b ⟺ name_a < name_b (alphabetically)
    all_var_names.sort();
    let var_to_id: HashMap<&str, u16> = all_var_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i as u16))
        .collect();
    let id_to_var: Vec<&str> = all_var_names.iter().map(|s| s.as_str()).collect();

    // Phase 2: Build per-star sorted variant ID arrays
    let mut dhap: HashMap<String, String> = HashMap::new();
    let mut dstar: HashMap<String, String> = HashMap::new();
    let mut star_var_ids: HashMap<String, Vec<u16>> = HashMap::new();

    for &(star_id, ref variant_list) in &star_entries {
        let mut var_ids: Vec<u16> = variant_list
            .iter()
            .map(|v| var_to_id[v])
            .collect();
        var_ids.sort_unstable();

        let var_key = merge_to_key(&[&var_ids], &id_to_var);
        dhap.entry(var_key.clone())
            .or_insert_with(|| star_id.to_string());
        dstar
            .entry(star_id.to_string())
            .or_insert_with(|| var_key);
        star_var_ids
            .entry(star_id.to_string())
            .or_insert(var_ids);
    }

    let mut star_keys: Vec<String> = star_var_ids.keys().cloned().collect();
    star_keys.sort();

    // Phase 3: Build combination tables using fast integer merging

    // CN=2: iterate ordered pairs (s1 <= s2) since result is symmetric
    let mut dhap2: HashMap<String, Vec<String>> = HashMap::new();
    let n = star_keys.len();
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in i..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            let key = merge_to_key(&[ids1, ids2], &id_to_var);
            // star_keys is sorted, so s1 <= s2 already — no need to sort
            let ss = format!("{}_{}", &star_keys[i], &star_keys[j]);
            insert_hap(&key, &ss, &mut dhap2);
        }
    }

    // CN=3: iterate ordered triples (s1 <= s2 <= s3)
    let mut dhap3: HashMap<String, Vec<String>> = HashMap::new();
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in i..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            for k in j..n {
                let ids3 = &star_var_ids[&star_keys[k]];
                let key = merge_to_key(&[ids1, ids2, ids3], &id_to_var);
                let ss = format!("{}_{}_{}", &star_keys[i], &star_keys[j], &star_keys[k]);
                insert_hap(&key, &ss, &mut dhap3);
            }
        }
    }

    // CN=3 pair (s1 + s2×2): iterate all pairs since s1 and s2 play different roles
    let mut dhap3pair: HashMap<String, Vec<String>> = HashMap::new();
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in 0..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            let key = merge_to_key(&[ids1, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap3pair);
        }
    }

    // CN=4 pair: s1×2+s2×2 (ordered i<=j) and s1+s2×3 (all pairs)
    let mut dhap4pair: HashMap<String, Vec<String>> = HashMap::new();
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in i..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            // s1×2 + s2×2
            let key = merge_to_key(&[ids1, ids1, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[i], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap4pair);
        }
        for j in 0..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            // s1 + s2×3
            let key = merge_to_key(&[ids1, ids2, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap4pair);
        }
    }

    // CN=5: s1×2+s2×3 (all pairs) and s1+s2×4 (all pairs)
    let mut dhap5pair: HashMap<String, Vec<String>> = HashMap::new();
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in 0..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            let key = merge_to_key(&[ids1, ids1, ids2, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap5pair);

            let key = merge_to_key(&[ids1, ids2, ids2, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap5pair);
        }
    }

    // CN=6: s1×2+s2×4, s1+s2×5, s1×3+s2×3
    let mut dhap6pair: HashMap<String, Vec<String>> = HashMap::new();
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in 0..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            let key = merge_to_key(&[ids1, ids1, ids2, ids2, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap6pair);

            let key = merge_to_key(&[ids1, ids2, ids2, ids2, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap6pair);
        }
    }
    for i in 0..n {
        let ids1 = &star_var_ids[&star_keys[i]];
        for j in i..n {
            let ids2 = &star_var_ids[&star_keys[j]];
            // s1×3 + s2×3 (ordered i<=j)
            let key = merge_to_key(&[ids1, ids1, ids1, ids2, ids2, ids2], &id_to_var);
            let ss = make_star_set(&[&star_keys[i], &star_keys[i], &star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j]]);
            insert_hap(&key, &ss, &mut dhap6pair);
        }
    }

    // exon9hybx2
    let mut dhap_exon9_x2: HashMap<String, Vec<String>> = HashMap::new();
    for &e1 in EXON9GC_ALLELES {
        for &e2 in EXON9GC_ALLELES {
            if let (Some(ids_e1), Some(ids_e2)) = (star_var_ids.get(e1), star_var_ids.get(e2)) {
                for s3 in &star_keys {
                    let ids3 = &star_var_ids[s3];
                    for s4 in &star_keys {
                        let ids4 = &star_var_ids[s4];
                        let key = merge_to_key(&[ids_e1, ids_e2, ids3, ids4], &id_to_var);
                        let ss = make_star_set(&[e1, e2, s3, s4]);
                        insert_hap(&key, &ss, &mut dhap_exon9_x2);
                    }
                }
            }
        }
    }

    // exon9hybx3
    let mut dhap_exon9_x3: HashMap<String, Vec<String>> = HashMap::new();
    for &e1 in EXON9GC_ALLELES {
        for &e2 in EXON9GC_ALLELES {
            for &e3 in EXON9GC_ALLELES {
                if let (Some(ids_e1), Some(ids_e2), Some(ids_e3)) =
                    (star_var_ids.get(e1), star_var_ids.get(e2), star_var_ids.get(e3))
                {
                    for s4 in &star_keys {
                        let ids4 = &star_var_ids[s4];
                        for s5 in &star_keys {
                            let ids5 = &star_var_ids[s5];
                            let key = merge_to_key(&[ids_e1, ids_e2, ids_e3, ids4, ids5], &id_to_var);
                            let ss = make_star_set(&[e1, e2, e3, s4, s5]);
                            insert_hap(&key, &ss, &mut dhap_exon9_x3);
                        }
                    }
                }
            }
        }
    }

    // exon9hybx4
    let mut dhap_exon9_x4: HashMap<String, Vec<String>> = HashMap::new();
    if let (Some(ids_10), Some(ids_36)) = (star_var_ids.get("*10"), star_var_ids.get("*36")) {
        let key = merge_to_key(&[ids_10, ids_10, ids_36, ids_36, ids_36, ids_36], &id_to_var);
        let ss = make_star_set(&["*10", "*10", "*36", "*36", "*36", "*36"]);
        insert_hap(&key, &ss, &mut dhap_exon9_x4);
    }

    // dup_exon9hyb
    let mut dhap_dup_exon9: HashMap<String, Vec<String>> = HashMap::new();
    let pair_alleles = exon9gc_pair_alleles();
    for (&e1, &e2) in &pair_alleles {
        if let (Some(ids_e1), Some(ids_e2)) = (star_var_ids.get(e1), star_var_ids.get(e2)) {
            for s3 in &star_keys {
                let ids3 = &star_var_ids[s3];
                for s4 in &star_keys {
                    let ids4 = &star_var_ids[s4];
                    let key = merge_to_key(&[ids_e1, ids_e2, ids3, ids4], &id_to_var);
                    let ss = make_star_set(&[e1, e2, s3, s4]);
                    insert_hap(&key, &ss, &mut dhap_dup_exon9);
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

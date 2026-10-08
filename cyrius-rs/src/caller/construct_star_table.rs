use crate::types::StarCombinations;
use std::collections::HashMap;
use std::thread::{Scope, ScopedJoinHandle};

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

/// A combination table: variant key -> the star-allele multisets that produce it.
pub type HapTable = HashMap<String, Vec<String>>;

/// Which of [`StarCombinations`]'s ten combination tables a CNV tag needs. Exactly
/// one is consulted per run, so this is also the unit of work for [`get_hap_table`]:
/// the other nine cost ~3 s and ~1 GB for nothing (finding #41).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableKind {
    Dhap2, Dhap3, Dhap3Pair, Dhap4Pair, Dhap5Pair,
    Dhap6Pair, Exon9X2, Exon9X3, Exon9X4, DupExon9,
}

impl TableKind {
    pub const ALL: [Self; 10] = [
        Self::Dhap2, Self::Dhap3, Self::Dhap3Pair, Self::Dhap4Pair, Self::Dhap5Pair,
        Self::Dhap6Pair, Self::Exon9X2, Self::Exon9X3, Self::Exon9X4, Self::DupExon9,
    ];

    /// The field this kind names: the one place field identity is decided.
    pub fn table(self, c: &StarCombinations) -> &HapTable {
        match self {
            TableKind::Dhap2 => &c.dhap2,
            TableKind::Dhap3 => &c.dhap3,
            TableKind::Dhap3Pair => &c.dhap3pair,
            TableKind::Dhap4Pair => &c.dhap4pair,
            TableKind::Dhap5Pair => &c.dhap5pair,
            TableKind::Dhap6Pair => &c.dhap6pair,
            TableKind::Exon9X2 => &c.dhap_exon9_x2,
            TableKind::Exon9X3 => &c.dhap_exon9_x3,
            TableKind::Exon9X4 => &c.dhap_exon9_x4,
            TableKind::DupExon9 => &c.dhap_dup_exon9,
        }
    }
}

/// The combination table a CNV/hybrid group matches against, or `None` for the
/// groups that match against `dhap` (CN=1) or do not match at all. The single copy
/// of this mapping: `match_star_allele::get_dic` resolves it to a reference and
/// [`get_hap_table`] reads it to decide what to build, so the two cannot drift --
/// which is what finding #26 was.
pub fn combination_table_for(cnvcall: &str) -> Option<TableKind> {
    match cnvcall {
        "star5" | "star13" | "star5_star68" | "star13_star68" => None, // Uses dhap (CN=1)
        "cn2" | "star13intron1" | "dup_star13" | "exon9hyb_star5" => Some(TableKind::Dhap2),
        s if s.starts_with("star68") && !s.contains("exon9hyb") && !s.contains("dup") => Some(TableKind::Dhap2),
        "exon9hyb" | "exon9hyb_star68" => Some(TableKind::Dhap3),
        "cn3" | "dup_star13intron1" | "dup_star68" => Some(TableKind::Dhap3Pair),
        "cn4" => Some(TableKind::Dhap4Pair),
        "cn5" => Some(TableKind::Dhap5Pair),
        "cn6" => Some(TableKind::Dhap6Pair),
        "exon9hyb_exon9hyb" => Some(TableKind::Exon9X2),
        "exon9hyb_exon9hyb_exon9hyb" => Some(TableKind::Exon9X3),
        "exon9hyb_exon9hyb_exon9hyb_exon9hyb" => Some(TableKind::Exon9X4),
        "dup_exon9hyb" => Some(TableKind::DupExon9),
        _ => None,
    }
}

/// Spawn a combination-table builder only if its table was asked for.
fn spawn_if<'s, 'e, F: FnOnce() -> HapTable + Send + 's>(
    s: &'s Scope<'s, 'e>, build: bool, f: F,
) -> Option<ScopedJoinHandle<'s, HapTable>> {
    build.then(|| s.spawn(f))
}

/// An unspawned builder yields the empty table.
fn opt_join(handle: Option<ScopedJoinHandle<'_, HapTable>>) -> HapTable {
    handle.map_or_else(HapTable::new, |h| h.join().unwrap())
}

/// Build the variant tables for the star allele definition file. `dhap` and `dstar`
/// are always built (cheap, and several callers read them); of the ten combination
/// tables only those in `wanted` are built, one scoped thread each, rest empty.
pub fn get_hap_table(star_table_content: &str, wanted: &[TableKind]) -> StarCombinations {
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

    // Phase 3: Build the wanted combination tables, one scoped thread each
    let n = star_keys.len();
    let want = |kind: TableKind| wanted.contains(&kind);

    let (dhap2, dhap3, dhap3pair, dhap4pair, dhap5pair, dhap6pair,
     dhap_exon9_x2, dhap_exon9_x3, dhap_exon9_x4, dhap_dup_exon9) =
        std::thread::scope(|s| {
            let h_dhap2 = spawn_if(s, want(TableKind::Dhap2), || {
                let mut dhap2: HashMap<String, Vec<String>> = HashMap::new();
                for i in 0..n {
                    let ids1 = &star_var_ids[&star_keys[i]];
                    for j in i..n {
                        let ids2 = &star_var_ids[&star_keys[j]];
                        let key = merge_to_key(&[ids1, ids2], &id_to_var);
                        let ss = format!("{}_{}", &star_keys[i], &star_keys[j]);
                        insert_hap(&key, &ss, &mut dhap2);
                    }
                }
                dhap2
            });

            let h_dhap3 = spawn_if(s, want(TableKind::Dhap3), || {
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
                dhap3
            });

            let h_dhap3pair = spawn_if(s, want(TableKind::Dhap3Pair), || {
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
                dhap3pair
            });

            let h_dhap4pair = spawn_if(s, want(TableKind::Dhap4Pair), || {
                let mut dhap4pair: HashMap<String, Vec<String>> = HashMap::new();
                for i in 0..n {
                    let ids1 = &star_var_ids[&star_keys[i]];
                    for j in i..n {
                        let ids2 = &star_var_ids[&star_keys[j]];
                        let key = merge_to_key(&[ids1, ids1, ids2, ids2], &id_to_var);
                        let ss = make_star_set(&[&star_keys[i], &star_keys[i], &star_keys[j], &star_keys[j]]);
                        insert_hap(&key, &ss, &mut dhap4pair);
                    }
                    for j in 0..n {
                        let ids2 = &star_var_ids[&star_keys[j]];
                        let key = merge_to_key(&[ids1, ids2, ids2, ids2], &id_to_var);
                        let ss = make_star_set(&[&star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j]]);
                        insert_hap(&key, &ss, &mut dhap4pair);
                    }
                }
                dhap4pair
            });

            let h_dhap5pair = spawn_if(s, want(TableKind::Dhap5Pair), || {
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
                dhap5pair
            });

            let h_dhap6pair = spawn_if(s, want(TableKind::Dhap6Pair), || {
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
                        let key = merge_to_key(&[ids1, ids1, ids1, ids2, ids2, ids2], &id_to_var);
                        let ss = make_star_set(&[&star_keys[i], &star_keys[i], &star_keys[i], &star_keys[j], &star_keys[j], &star_keys[j]]);
                        insert_hap(&key, &ss, &mut dhap6pair);
                    }
                }
                dhap6pair
            });

            let h_exon9_x2 = spawn_if(s, want(TableKind::Exon9X2), || {
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
                dhap_exon9_x2
            });

            let h_exon9_x3 = spawn_if(s, want(TableKind::Exon9X3), || {
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
                dhap_exon9_x3
            });

            let h_exon9_x4 = spawn_if(s, want(TableKind::Exon9X4), || {
                let mut dhap_exon9_x4: HashMap<String, Vec<String>> = HashMap::new();
                if let (Some(ids_10), Some(ids_36)) = (star_var_ids.get("*10"), star_var_ids.get("*36")) {
                    let key = merge_to_key(&[ids_10, ids_10, ids_36, ids_36, ids_36, ids_36], &id_to_var);
                    let ss = make_star_set(&["*10", "*10", "*36", "*36", "*36", "*36"]);
                    insert_hap(&key, &ss, &mut dhap_exon9_x4);
                }
                dhap_exon9_x4
            });

            let h_dup_exon9 = spawn_if(s, want(TableKind::DupExon9), || {
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
                dhap_dup_exon9
            });

            (opt_join(h_dhap2), opt_join(h_dhap3), opt_join(h_dhap3pair),
             opt_join(h_dhap4pair), opt_join(h_dhap5pair), opt_join(h_dhap6pair),
             opt_join(h_exon9_x2), opt_join(h_exon9_x3), opt_join(h_exon9_x4),
             opt_join(h_dup_exon9))
        });

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
        dstar,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::caller::call_variants::CNV_ACCEPTED;
    use crate::data::STAR_TABLE;

    /// The property the lazy build rests on (finding #41): for every accepted CNV
    /// tag, a build of just the table that tag selects must put byte-identical
    /// content in that field, leave the other nine empty, and still carry the two
    /// cheap always-built tables.
    ///
    /// The nine-empty half is not pedantry: `get_dic` reads whichever field the tag
    /// selects, so a kind written into the wrong field would hand the matcher an
    /// empty map and silently drop every candidate.
    #[test]
    fn lazy_build_matches_eager_build_for_every_accepted_tag() {
        let eager = get_hap_table(STAR_TABLE, &TableKind::ALL);
        for tag in CNV_ACCEPTED {
            let kind = combination_table_for(tag);
            let lazy = get_hap_table(STAR_TABLE, kind.as_slice());
            assert_eq!(lazy.dhap, eager.dhap, "dhap differs for tag {}", tag);
            assert_eq!(lazy.dstar, eager.dstar, "dstar differs for tag {}", tag);
            for other in TableKind::ALL {
                if Some(other) == kind {
                    // `dhap_dup_exon9` is the one table whose per-key candidate
                    // *order* is not reproducible even between two eager builds:
                    // `exon9gc_pair_alleles` returns a `HashMap`, and each
                    // `HashMap` instance gets its own `RandomState`, so the two
                    // pairs are visited in an arbitrary order and `insert_hap`
                    // appends in that order. Compare its content, not its order;
                    // the other nine are compared exactly.
                    let mut got = other.table(&lazy).clone();
                    let mut expected = other.table(&eager).clone();
                    if other == TableKind::DupExon9 {
                        got.values_mut().for_each(|v| v.sort());
                        expected.values_mut().for_each(|v| v.sort());
                    }
                    assert_eq!(
                        got, expected,
                        "{:?} differs from the eager build for tag {}",
                        other, tag
                    );
                    assert!(
                        !other.table(&lazy).is_empty(),
                        "{:?} is the table tag {} selects, yet it was not built",
                        other,
                        tag
                    );
                } else {
                    assert!(
                        other.table(&lazy).is_empty(),
                        "tag {} selects {:?}, so {:?} should not have been built",
                        tag,
                        kind,
                        other
                    );
                }
            }
        }
    }

    /// Every kind must be selected by some accepted tag. A kind nothing selects
    /// would be a table built by nothing and checked by nothing above.
    #[test]
    fn every_table_kind_is_selected_by_some_accepted_tag() {
        for kind in TableKind::ALL {
            assert!(
                CNV_ACCEPTED
                    .iter()
                    .any(|tag| combination_table_for(tag) == Some(kind)),
                "no accepted CNV tag selects {:?}",
                kind
            );
        }
    }
}

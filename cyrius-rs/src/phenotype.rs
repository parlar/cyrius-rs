use crate::types::{FrequencyRow, HaplotypeDetail, HaplotypeFunctionality, Prediction};
use std::collections::HashMap;

/// Load CYP2D6 haplotype functionality from embedded data.
pub fn load_haplotype_functionality(content: &str) -> HashMap<String, HaplotypeFunctionality> {
    let mut result = HashMap::new();
    let mut lines = content.lines();
    let _header = lines.next(); // Skip header

    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let haplotype = fields[0].to_string();
        let activity = fields[1].to_string();
        let function = fields[2].to_string();
        let evidence_strength = fields[3].to_string();
        let evidence_summary = fields[4].to_string();
        result.insert(
            haplotype,
            HaplotypeFunctionality {
                activity,
                function,
                evidence_strength,
                evidence_summary,
            },
        );
    }
    result
}

/// Format an activity value the way Python's str(float) does.
fn format_activity(value: f64) -> String {
    if value == value.floor() {
        format!("{:.1}", value)
    } else {
        format!("{}", value)
    }
}

/// Look up a haplotype's functionality. When there is no exact row:
/// a sub-allele uses its main allele ("*4.013" -> "*4"); otherwise "*AxN"
/// gets N times the single-copy activity of "*A". The table's "*Ax≥3" rows
/// are lower bounds for an unknown count, so they only give the function
/// label for N >= 3.
fn lookup_functionality(
    name: &str,
    table: &HashMap<String, HaplotypeFunctionality>,
) -> Option<HaplotypeFunctionality> {
    if let Some(info) = table.get(name) {
        return Some(info.clone());
    }
    let (base, copies) = match name.rsplit_once('x') {
        Some((b, n)) if !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()) => {
            (b, n.parse::<u32>().ok())
        }
        _ => (name, None),
    };
    let main = base.split('.').next().unwrap_or(base);
    if main != base {
        let alt = match copies {
            Some(n) => format!("{}x{}", main, n),
            None => main.to_string(),
        };
        return lookup_functionality(&alt, table);
    }
    let n = copies?;
    let single = table.get(base)?;
    let activity: f64 = single.activity.parse().ok()?;
    // The copy count is exact, so the activity is n times one copy. A
    // "*Ax≥3" row is a lower bound for an unknown count of 3 or more: it
    // supplies only the function label, never the activity.
    let function = match table.get(&format!("{}x≥3", base)) {
        Some(bucket) if n >= 3 => bucket.function.clone(),
        _ => single.function.clone(),
    };
    Some(HaplotypeFunctionality {
        activity: format_activity(activity * n as f64),
        function,
        evidence_strength: "n/a".to_string(),
        evidence_summary: format!(
            "No table entry for {}; activity is {} x the {} value.",
            name, n, base
        ),
    })
}

/// Match haplotype with activity value, functionality, evidence.
pub fn compare_genotypes(
    query: &str,
    haplotype_functionality: &HashMap<String, HaplotypeFunctionality>,
) -> (Vec<HaplotypeDetail>, bool) {
    let mut haplotypes = Vec::new();
    for haplotype_pair in query.split('/') {
        for subhaplotype in haplotype_pair.split('+') {
            haplotypes.push(subhaplotype.trim().to_string());
        }
    }

    let mut details = Vec::new();
    let mut has_na = false;

    for haplotype in &haplotypes {
        if let Some(info) = lookup_functionality(haplotype, haplotype_functionality) {
            if info.activity == "n/a" {
                has_na = true;
            }
            details.push(HaplotypeDetail {
                haplotype: haplotype.clone(),
                activity: info.activity.clone(),
                function: info.function.clone(),
                evidence_strength: info.evidence_strength.clone(),
                evidence_summary: info.evidence_summary.clone(),
            });
        } else {
            has_na = true;
            details.push(HaplotypeDetail {
                haplotype: haplotype.clone(),
                activity: "n/a".to_string(),
                function: "n/a".to_string(),
                evidence_strength: "n/a".to_string(),
                evidence_summary: "n/a".to_string(),
            });
        }
    }

    (details, has_na)
}

/// Calculate total activity score from haplotype details.
pub fn calculate_total_activity_score(haplotype_details: &[HaplotypeDetail]) -> String {
    let mut total = 0.0f64;
    for detail in haplotype_details {
        if detail.activity == "n/a" {
            return "n/a".to_string();
        }
        let cleaned = detail.activity.replace("≥", "");
        match cleaned.parse::<f64>() {
            Ok(v) => total += v,
            Err(_) => return "n/a".to_string(),
        }
    }
    // Match Python str(float) behavior: always show at least one decimal place
    if total == total.floor() {
        format!("{:.1}", total)
    } else {
        format!("{}", total)
    }
}

/// Determine phenotype based on activity score.
pub fn determine_phenotype(activity_score: &str) -> String {
    if activity_score == "n/a" {
        return "Indeterminate".to_string();
    }
    match activity_score.parse::<f64>() {
        Ok(score) if score.is_finite() => {
            // Contiguous thresholds (no gaps): identical classification to the
            // CPIC 0.25-increment bands for valid scores, but fractional values
            // that fall between increments still resolve sensibly.
            if score <= 0.0 {
                "Poor Metabolizer".to_string()
            } else if score <= 1.0 {
                "Intermediate Metabolizer".to_string()
            } else if score <= 2.25 {
                "Normal Metabolizer".to_string()
            } else {
                "Ultrarapid Metabolizer".to_string()
            }
        }
        _ => "Unknown".to_string(),
    }
}

/// Sort genotype so smaller integer will be the first element.
pub fn sort_genotype(genotype: Option<&str>) -> Option<String> {
    let genotype = genotype?;
    if genotype.is_empty() {
        return Some(genotype.to_string());
    }

    let parts: Vec<&str> = genotype.split('/').collect();
    if parts.len() != 2 {
        return Some(genotype.to_string());
    }

    let parse_star_num = |s: &str| -> i64 {
        let after_star = s.split('*').nth(1).unwrap_or("");
        let num_str = if after_star.contains('+') {
            after_star.split('+').next().unwrap_or("")
        } else {
            after_star.split('x').next().unwrap_or("")
        };
        num_str.parse::<i64>().unwrap_or(i64::MAX)
    };

    let a1_int = parse_star_num(parts[0]);
    let a2_int = parse_star_num(parts[1]);

    if a1_int > a2_int {
        Some(format!("{}/{}", parts[1], parts[0]))
    } else {
        Some(genotype.to_string())
    }
}

/// Compare the predicted effect of ';'-separated alternative diplotypes.
/// Returns (clinical_equivalence, long_read_recommended). Any alternative
/// whose activity is unknown makes the effect unknown, which also warrants
/// long-read follow-up.
pub fn clinical_equivalence(predictions: &[Prediction]) -> (String, bool) {
    if predictions.iter().any(|p| p.total_activity == "n/a") {
        return ("unknown_effect".to_string(), true);
    }
    let all_same = predictions.windows(2).all(|w| {
        (w[0].total_activity.as_str(), w[0].predicted_phenotype.as_str())
            == (w[1].total_activity.as_str(), w[1].predicted_phenotype.as_str())
    });
    if all_same {
        ("same_effect".to_string(), false)
    } else {
        ("different_effect".to_string(), true)
    }
}

/// Apply sort_genotype to each ';'-separated diplotype.
pub fn sort_diplotypes(genotype: &str) -> String {
    genotype
        .split(';')
        .map(|d| sort_genotype(Some(d.trim())).unwrap_or_default())
        .collect::<Vec<_>>()
        .join(";")
}

/// Match phenotype for a genotype string (may contain multiple diplotypes separated by ;).
pub fn match_phenotype(
    genotype: Option<&str>,
    haplotype_functionality: &HashMap<String, HaplotypeFunctionality>,
) -> Vec<Prediction> {
    let genotype = match genotype {
        Some(g) if !g.is_empty() && g != "None" => g,
        _ => {
            return vec![Prediction {
                total_activity: "n/a".to_string(),
                predicted_phenotype: "n/a".to_string(),
                haplotype_details: vec![],
            }]
        }
    };

    let mut predictions = Vec::new();
    for diplotype in genotype.split(';') {
        let sorted = sort_genotype(Some(diplotype.trim())).unwrap_or_default();
        let (details, has_na) = compare_genotypes(&sorted, haplotype_functionality);
        let total_activity = if has_na {
            "n/a".to_string()
        } else {
            calculate_total_activity_score(&details)
        };
        let predicted_phenotype = determine_phenotype(&total_activity);
        predictions.push(Prediction {
            total_activity,
            predicted_phenotype,
            haplotype_details: details,
        });
    }

    if predictions.is_empty() {
        vec![Prediction {
            total_activity: "n/a".to_string(),
            predicted_phenotype: "n/a".to_string(),
            haplotype_details: vec![],
        }]
    } else {
        predictions
    }
}

/// Get population frequency percentages for a haplotype.
pub fn get_percentages(content: &str, name: &str) -> Option<Vec<(String, String)>> {
    let mut lines = content.lines();
    let header_line = lines.next()?;
    let headers: Vec<&str> = header_line.split('\t').collect();

    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        if fields[0] == name {
            let mut result = Vec::new();
            for (i, &header) in headers.iter().enumerate().skip(1) {
                if i < fields.len() {
                    let val: f64 = fields[i].parse().unwrap_or(f64::NAN);
                    if val.is_nan() {
                        result.push((header.to_string(), String::new()));
                    } else {
                        let pct = (val * 100.0 * 1000.0).round() / 1000.0;
                        // Match Python str(float) behavior: always show at least one decimal
                        let pct_str = if pct == pct.floor() {
                            format!("{:.1}%", pct)
                        } else {
                            format!("{}%", pct)
                        };
                        result.push((header.to_string(), pct_str));
                    }
                }
            }
            return Some(result);
        }
    }
    None
}

/// Build diplotype frequency table.
pub fn diplotype_frequencies(
    haplotype_content: &str,
    diplotype_content: &str,
    diplotype_string: &str,
) -> Result<Vec<FrequencyRow>, String> {
    let mut all_columns: Vec<String> = Vec::new();
    let mut column_data: indexmap::IndexMap<String, indexmap::IndexMap<String, String>> =
        indexmap::IndexMap::new();

    for diplotype in diplotype_string.split(';') {
        let dip_pcts = get_percentages(diplotype_content, diplotype)
            .ok_or_else(|| format!("Diplotype {} not found in the file.", diplotype))?;
        all_columns.push(diplotype.to_string());
        let mut col = indexmap::IndexMap::new();
        for (group, val) in &dip_pcts {
            col.insert(group.clone(), val.clone());
        }
        column_data.insert(diplotype.to_string(), col);

        for haplotype in diplotype.split('/') {
            if column_data.contains_key(haplotype) {
                continue;
            }
            let hap_pcts = get_percentages(haplotype_content, haplotype)
                .ok_or_else(|| format!("Haplotype {} not found in the file.", haplotype))?;
            all_columns.push(haplotype.to_string());
            let mut col = indexmap::IndexMap::new();
            for (group, val) in &hap_pcts {
                col.insert(group.clone(), val.clone());
            }
            column_data.insert(haplotype.to_string(), col);
        }
    }

    // Get biogeographic groups from first column
    let first_col = column_data
        .values()
        .next()
        .ok_or_else(|| "No diplotype/haplotype frequency data available".to_string())?;
    let groups: Vec<String> = first_col.keys().cloned().collect();

    let mut rows = Vec::new();
    for group in &groups {
        let mut values = indexmap::IndexMap::new();
        for col_name in &all_columns {
            let val = column_data[col_name]
                .get(group)
                .cloned()
                .unwrap_or_default();
            values.insert(col_name.clone(), val);
        }
        rows.push(FrequencyRow {
            biogeographic_group: group.clone(),
            values,
        });
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> HashMap<String, HaplotypeFunctionality> {
        load_haplotype_functionality(crate::data::HAPLOTYPE_FUNC)
    }

    fn activity_of(query: &str, hap: &str) -> String {
        let (details, _) = compare_genotypes(query, &table());
        details.into_iter().find(|d| d.haplotype == hap).unwrap().activity
    }

    #[test]
    fn duplication_without_its_own_row_multiplies_single_copy_activity() {
        assert_eq!(activity_of("*1/*10x3", "*10x3"), "0.75");
    }

    #[test]
    fn an_exact_copy_count_multiplies_the_single_copy_activity() {
        // *41 is 0.25 per copy; the *41x3 row is an exact 3-copy row, not a
        // ">= 3" bucket, so 4 and 5 copies must not reuse its 0.75.
        assert_eq!(activity_of("*10/*41x4", "*41x4"), "1.0");
        assert_eq!(activity_of("*5/*41x5", "*41x5"), "1.25");
        // *2x≥3 is a lower bound ("≥3.0"); an exact count of 4 is 4.0.
        assert_eq!(activity_of("*5/*2x4", "*2x4"), "4.0");
    }

    #[test]
    fn an_exact_table_row_still_wins() {
        assert_eq!(activity_of("*1/*41x3", "*41x3"), "0.75");
        assert_eq!(activity_of("*1/*2x2", "*2x2"), "2.0");
    }

    #[test]
    fn exact_multiplication_sets_the_cpic_phenotype() {
        let t = table();
        let p = match_phenotype(Some("*41x4/*10"), &t);
        assert_eq!(p[0].total_activity, "1.25");
        assert_eq!(p[0].predicted_phenotype, "Normal Metabolizer");
        let p = match_phenotype(Some("*41x5/*5"), &t);
        assert_eq!(p[0].total_activity, "1.25");
        assert_eq!(p[0].predicted_phenotype, "Normal Metabolizer");
    }

    #[test]
    fn three_or_more_normal_copies_keep_the_increased_function_label() {
        let (details, _) = compare_genotypes("*1x4/*5", &table());
        let d = details.into_iter().find(|d| d.haplotype == "*1x4").unwrap();
        assert_eq!(d.activity, "4.0");
        assert_eq!(d.function, "Increased function");
    }

    #[test]
    fn suballele_falls_back_to_its_main_allele() {
        assert_eq!(activity_of("*4.013/*5", "*4.013"), "0.0");
    }

    fn pred(activity: &str, phenotype: &str) -> Prediction {
        Prediction {
            total_activity: activity.to_string(),
            predicted_phenotype: phenotype.to_string(),
            haplotype_details: vec![],
        }
    }

    #[test]
    fn unknown_alternatives_recommend_long_read() {
        let p = [pred("n/a", "Indeterminate"), pred("n/a", "Indeterminate")];
        assert_eq!(clinical_equivalence(&p), ("unknown_effect".to_string(), true));
        let same = [pred("1.0", "Intermediate Metabolizer"), pred("1.0", "Intermediate Metabolizer")];
        assert_eq!(clinical_equivalence(&same), ("same_effect".to_string(), false));
        let diff = [pred("0.25", "Intermediate Metabolizer"), pred("0.5", "Intermediate Metabolizer")];
        assert_eq!(clinical_equivalence(&diff), ("different_effect".to_string(), true));
    }

    #[test]
    fn each_diplotype_is_sorted() {
        assert_eq!(sort_diplotypes("*2/*1;*10/*4"), "*1/*2;*4/*10");
    }
}

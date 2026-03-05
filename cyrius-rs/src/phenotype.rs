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
        let haplotype = fields[0].replace("≥", "");
        let activity = fields[1].replace("≥", "");
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
        if let Some(info) = haplotype_functionality.get(haplotype) {
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
        Ok(score) => {
            if score == 0.0 {
                "Poor Metabolizer".to_string()
            } else if (0.25..=1.0).contains(&score) {
                "Intermediate Metabolizer".to_string()
            } else if (1.25..=2.25).contains(&score) {
                "Normal Metabolizer".to_string()
            } else if score > 2.25 {
                "Ultrarapid Metabolizer".to_string()
            } else {
                "Unknown".to_string()
            }
        }
        Err(_) => "Unknown".to_string(),
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
    let first_col = column_data.values().next().unwrap();
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

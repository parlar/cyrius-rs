use crate::depth_calling::snp_count::passing_read;
use crate::depth_calling::utilities::fetch_region;
use crate::types::SnpLookup;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

/// Get haplotypes for each read at a set of target positions (both regions).
pub fn get_haplotypes_from_bam(
    reader: &mut bam::IndexedReader,
    base_db: &SnpLookup,
    target_positions: &[usize],
) -> Result<HashMap<String, String>, String> {
    let dread = get_bases_per_read(reader, base_db, target_positions, None, 0)?;
    let (base1, base2) = get_base1_base2(base_db, target_positions);
    Ok(get_hap_counts(&dread, &base1, &base2, target_positions))
}

/// Get haplotypes for each read at a set of target positions (single region only).
pub fn get_haplotypes_from_bam_single_region(
    reader: &mut bam::IndexedReader,
    base_db: &SnpLookup,
    target_positions: &[usize],
) -> Result<HashMap<String, String>, String> {
    let dread = get_bases_per_read(reader, base_db, target_positions, Some(0), 10)?;
    let (base1, base2) = get_base1_base2(base_db, target_positions);
    Ok(get_hap_counts(&dread, &base1, &base2, target_positions))
}

/// Get bases per read at target positions from pileup data.
/// region: None = both, Some(0) = dsnp1 only, Some(1) = dsnp2 only.
fn get_bases_per_read(
    reader: &mut bam::IndexedReader,
    base_db: &SnpLookup,
    target_positions: &[usize],
    region: Option<usize>,
    min_mapq: u8,
) -> Result<HashMap<String, HashMap<usize, Option<String>>>, String> {
    let mut dread: HashMap<String, HashMap<usize, Option<String>>> = HashMap::new();
    let nchr = &base_db.nchr;
    let dindex = &base_db.dindex;

    let dsnps: Vec<&indexmap::IndexMap<String, String>> = match region {
        Some(0) => vec![&base_db.dsnp1],
        Some(1) => vec![&base_db.dsnp2],
        _ => vec![&base_db.dsnp1, &base_db.dsnp2],
    };

    let tid = match reader.header().tid(nchr.as_bytes()) {
        Some(t) => t,
        None => return Ok(dread),
    };

    for dsnp in &dsnps {
        for (snp_position_ori, alleles) in dsnp.iter() {
            let dsnp_index = match dindex.get(snp_position_ori) {
                Some(&v) => v,
                None => continue,
            };
            let snp_position: i64 = match snp_position_ori.split('_').next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };

            if !target_positions.contains(&dsnp_index) {
                continue;
            }

            // A failed fetch is not "no reads here": the reads were never
            // looked at, so the caller must no-call (E1).
            fetch_region(reader, tid, snp_position - 1, snp_position, &format!("{nchr}:{snp_position}"))?;

            let mut pileups = reader.pileup();
            pileups.set_max_depth(i32::MAX as u32);

            for pileup_result in pileups {
                let pileup = match pileup_result {
                    Ok(p) => p,
                    Err(e) => return Err(format!("BAM pileup failed while reading haplotypes: {}", e)),
                };
                let site_position = pileup.pos() as i64 + 1;
                if site_position != snp_position {
                    continue;
                }

                let (reg1_allele, reg2_allele) = match alleles.split_once('_') {
                    Some(pair) => pair,
                    None => continue,
                };

                for alignment in pileup.alignments() {
                    let record = alignment.record();
                    let is_del = alignment.is_del();
                    let is_refskip = alignment.is_refskip();

                    if !passing_read(&record, is_del, is_refskip) || record.mapq() < min_mapq {
                        continue;
                    }

                    let qpos = match alignment.qpos() {
                        Some(p) => p,
                        None => continue,
                    };

                    let read_name = String::from_utf8_lossy(record.qname()).to_string();
                    let seq = record.seq().as_bytes();

                    // `min_len` is the length of the ONE slice extracted for
                    // this (read, site) -- and therefore the only length at
                    // which `get_hap_counts` can ever recognise an allele.
                    //
                    // Both regions' allele fields can hold several
                    // comma-separated alleles ('A,G'), so the length of the
                    // JOINED `reg1_allele`/`reg2_allele` string is not an
                    // allele length at all: it counts the commas and the extra
                    // alleles. Minimise over the comma-SPLIT alleles of BOTH
                    // regions instead -- the same split `get_hap_counts` uses
                    // -- so the extraction length is always a real allele
                    // length (review finding #7).
                    //
                    // What the raw string length did wrong. Let `l_min` be the
                    // minimum length over all comma-split alleles of both
                    // regions, and `raw = min(|reg1_allele|, |reg2_allele|)`.
                    // The two are equal IFF some region holds exactly one
                    // allele whose length is `l_min`; otherwise `raw > l_min`
                    // and the over-long slice either
                    //   * matches no allele at all, so the read scores 'x' and
                    //     `extract_hap` discards it, or
                    //   * coincides with the length of a DIFFERENT allele at
                    //     the same site, silently cross-assigning the read to
                    //     the wrong region.
                    // An earlier revision of this comment claimed the first
                    // outcome always held. It does not, and the counterexample
                    // is this commit's own fixture ('A,G' vs 'AC', in
                    // `characterization_per_allele_length_fixes_commas_but_not_unequal_lengths`
                    // below): there the raw bound of 2 scores the read '2', not
                    // 'x'.
                    //
                    // The case that makes the split version STRICTLY better:
                    // both regions comma'd with all alleles the same length,
                    // e.g. 'A,G' vs 'C,T'. `raw` is then 3 -- a length no
                    // allele has -- so EVERY read at the site scores 'x', the
                    // panel reads as "no reads at all", and every downstream
                    // gate defaults false; a silent, total loss of the site.
                    // Verified by executing this function in
                    // `both_regions_comma_separated_with_equal_length_alleles_are_still_scored`
                    // below: raw yields {carries_a: "x", carries_t: "x"} and an
                    // EMPTY `extract_hap`, split yields {carries_a: "1",
                    // carries_t: "2"}. The case where NEITHER value is right is
                    // alleles of unequal length -- the residual limitation at
                    // the end of this comment.
                    //
                    // Real data, all four shipped position DBs (285 rows):
                    // ZERO rows have both regions comma'd; `raw` and `l_min`
                    // differ on exactly 5 rows, all indel/repeat rows of
                    // `data/CYP2D6_target_variant_38.txt` (lines 47, 49, 67,
                    // 80, 81), none of which reaches `get_bases_per_read`; and
                    // on all 5 the `l_min` allele is not a prefix of any longer
                    // allele at its site, so both bounds fail closed there.
                    // `data/CYP2D6_haplotype_38.txt` -- the only DB that
                    // reaches this function -- has ZERO commas, so the split is
                    // a verified strict no-op on reachable data.
                    //
                    // NOT yet consistent with `snp_count`, contrary to the plan
                    // for this change ("per-allele lengths (split on ','),
                    // matching `snp_count`"). Only the comparison OPERATOR now
                    // matches: `snp_count`'s `qpos + max_allele_len >
                    // seq.len()` reject and the `end_pos > seq.len()` reject
                    // below are the same test. The AGGREGATION still differs --
                    // `snp_count` takes the MAX over all of a site's alleles
                    // and then matches each allele at its own length, while
                    // this takes the MIN -- so `snp_count` requires room for
                    // the LONGEST allele and this requires room only for the
                    // SHORTEST, and at any site whose alleles differ in length
                    // the two functions still admit different read sets. They
                    // agree today only because all EIGHT sites reachable from
                    // here have equal reg1/reg2 lengths:
                    // `data/CYP2D6_haplotype_38.txt` lines 4-9 are 1/1, line 2
                    // is 10/10, line 3 is 12/12. `min` is deliberately NOT
                    // changed to `max`: that would move read admission with no
                    // mandate, and the benchmark has no power to detect it.
                    //
                    // RESIDUAL LIMITATION, not fixed here: one fixed-length
                    // slice per (read, site) cannot discriminate alleles of
                    // unequal length the way `snp_count` does by matching each
                    // allele at its own length. Unreachable from here today --
                    // the haplotype DB's single unequal-length row (line 10,
                    // g.42130653 `GCAGG` vs `GCAAGG`) is consumed only by
                    // `call_var42130655ins_a` through
                    // `get_supporting_reads_single_region`.
                    let min_len = reg1_allele
                        .split(',')
                        .chain(reg2_allele.split(','))
                        .map(str::len)
                        .min()
                        .unwrap_or(0);

                    // Fail CLOSED on a zero-length slice, which the comma split
                    // above made reachable. `min_len == 0` means some allele
                    // field held an EMPTY comma component -- 'A,' splits to
                    // ["A", ""] -- i.e. one trailing-comma typo in a
                    // hand-maintained 10-line TSV, in a project whose sibling
                    // DBs already carry 22 comma'd rows (21 in
                    // `data/CYP2D6_target_variant_38.txt`, 1 in
                    // `data/CYP2D6_SNP_38.txt`).
                    //
                    // Without this guard the failure is NOT a dropped read, it
                    // is a FABRICATED match on every read: `hap` would be "",
                    // `get_hap_counts` splits `base1`/`base2` on ',' too, and
                    // `"" == ""` MATCHES -- so every covering read at the site
                    // scores '1' (or '2') and is counted by `extract_hap`.
                    // Pinned in
                    // `a_trailing_comma_must_not_fabricate_a_match_on_every_read`.
                    //
                    // This is strictly NEW exposure, not a restatement of the
                    // old one: the raw `min(|reg1|, |reg2|)` needed an empty
                    // RAW field, which `get_snp_position` cannot emit -- its
                    // fields come from `split_whitespace`, so none is ever
                    // empty. The comma split needs only an empty COMPONENT.
                    if min_len == 0 {
                        continue;
                    }

                    let end_pos = qpos + min_len;
                    // `>`, not `>=`: `seq[qpos..end_pos]` is a valid slice when
                    // `end_pos == seq.len()`, i.e. when the extracted bases run
                    // to the read's last stored base. The old `>=` dropped
                    // those reads from the haplotypes while `snp_count`'s bound
                    // (`qpos + max_allele_len > seq.len()`) still counted the
                    // very same reads in the D6/D7 depth ratios -- the two
                    // copies of this pileup loop had drifted apart here.
                    //
                    // UPSTREAM DIVERGENCE, disclosed: `bcyrius`'s
                    // `depth_calling/haplotype.py` reads
                    //     end_pos = start_pos + min(len(reg1_allele), len(reg2_allele))
                    //     if end_pos < len(read_seq):
                    // so the pre-change Rust was a FAITHFUL port of upstream on
                    // both points, and upstream's own haplotype.py /
                    // snp_count.py pair carries the same min-vs-max
                    // inconsistency described above. Haplotype counts out of
                    // this function will therefore differ from published Cyrius
                    // output.
                    if end_pos > seq.len() {
                        continue;
                    }

                    let hap = String::from_utf8_lossy(&seq[qpos..end_pos]).to_string();

                    let entry = dread.entry(read_name).or_insert_with(|| {
                        let mut m = HashMap::new();
                        for &pos in target_positions {
                            m.insert(pos, None);
                        }
                        m
                    });

                    let current = entry.get(&dsnp_index).cloned().flatten();
                    if current.is_some() && current.as_deref() != Some(&hap) {
                        // Conflicting base — set to None
                        entry.insert(dsnp_index, None);
                    } else {
                        entry.insert(dsnp_index, Some(hap));
                    }
                }
            }
        }
    }
    Ok(dread)
}

/// Translate bases into haplotype strings (1 = reg1 allele, 2 = reg2 allele, x = unknown).
fn get_hap_counts(
    dread: &HashMap<String, HashMap<usize, Option<String>>>,
    base1: &[Option<String>],
    base2: &[Option<String>],
    target_positions: &[usize],
) -> HashMap<String, String> {
    let mut dread_hap = HashMap::new();
    for (read_name, read_bases) in dread {
        let mut pos_list: Vec<char> = vec!['x'; target_positions.len()];
        for (i, &pos) in target_positions.iter().enumerate() {
            if let Some(base_opt) = read_bases.get(&pos) {
                if let Some(ref base) = base_opt {
                    if let Some(ref b1) = base1[i] {
                        for allele in b1.split(',') {
                            if base == &allele.to_uppercase() {
                                pos_list[i] = '1';
                            }
                        }
                    }
                    if let Some(ref b2) = base2[i] {
                        for allele in b2.split(',') {
                            if base == &allele.to_uppercase() {
                                pos_list[i] = '2';
                            }
                        }
                    }
                }
            }
        }
        dread_hap.insert(read_name.clone(), pos_list.iter().collect());
    }
    dread_hap
}

/// Get expected bases at target positions for each allele/paralog.
fn get_base1_base2(
    base_db: &SnpLookup,
    target_positions: &[usize],
) -> (Vec<Option<String>>, Vec<Option<String>>) {
    let n = base_db.dsnp1.len();
    let mut base1: Vec<Option<String>> = vec![None; n];
    let mut base2: Vec<Option<String>> = vec![None; n];

    for (pos, alleles) in &base_db.dsnp1 {
        let dsnp_index = match base_db.dindex.get(pos) {
            Some(&v) => v,
            None => continue,
        };
        if target_positions.contains(&dsnp_index) {
            let index: usize = match pos.split('_').nth(1).and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            if index >= n {
                continue;
            }
            let (allele1, allele2) = match alleles.split_once('_') {
                Some(pair) => pair,
                None => continue,
            };
            base1[index] = Some(allele1.to_string());
            base2[index] = Some(allele2.to_string());
        }
    }

    (
        target_positions.iter().map(|&i| base1[i].clone()).collect(),
        target_positions.iter().map(|&i| base2[i].clone()).collect(),
    )
}

/// Extract haplotypes at certain positions and count them.
pub fn extract_hap(
    dhaplotype: &HashMap<String, String>,
    positions_to_extract: &[usize],
) -> HashMap<String, Vec<u32>> {
    let mut hap_count: HashMap<String, Vec<u32>> = HashMap::new();
    for hap in dhaplotype.values() {
        let hap_chars: Vec<char> = hap.chars().collect();
        let hap_base: Vec<char> = positions_to_extract.iter().map(|&p| hap_chars[p]).collect();
        if hap_base.contains(&'x') {
            continue;
        }
        let key: String = hap_base.iter().collect();
        hap_count.entry(key).or_insert_with(Vec::new).push(1);
    }
    hap_count
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One fixture read: `(qname, 0-based start, [(qpos, bases)])`.
    type FixtureRead<'a> = (&'a str, i64, &'a [(usize, &'a str)]);

    /// Write a small indexed BAM of all-match reads.
    ///
    /// Each read is `(qname, 0-based start, [(qpos, bases)])`; the rest of the
    /// read is `N`, which matches no allele used below. Reads must be given in
    /// coordinate order so the BAI build succeeds.
    ///
    /// NOTE: `bam::Record::new()` leaves `flags = 4` (unmapped), and htslib's
    /// `bam_endpos` ignores the CIGAR for an unmapped record, so an unflagged
    /// record spans only `[pos, pos + 1)` and silently drops out of the region
    /// queries `get_bases_per_read` issues. The `set_flags(0)` below is
    /// load-bearing.
    fn write_hap_bam(
        path: &std::path::Path,
        read_len: usize,
        reads: &[FixtureRead],
    ) {
        let mut header = bam::header::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr22");
        sq.push_tag(b"LN", 50_818_468);
        header.push_record(&sq);

        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            let cigar =
                bam::record::CigarString(vec![bam::record::Cigar::Match(read_len as u32)]);
            let qual = vec![30u8; read_len];
            for (name, pos, bases) in reads {
                let mut seq = vec![b'N'; read_len];
                for (qpos, allele) in bases.iter() {
                    seq[*qpos..*qpos + allele.len()].copy_from_slice(allele.as_bytes());
                }
                let mut rec = bam::Record::new();
                rec.set(name.as_bytes(), Some(&cigar), &seq, &qual);
                rec.set_tid(0);
                rec.set_pos(*pos);
                rec.set_mapq(60);
                rec.set_flags(0); // mapped, primary, forward, not a duplicate
                writer.write(&rec).unwrap();
            }
        }
        bam::index::build(path, None::<&std::path::Path>, bam::index::Type::Bai, 1).unwrap();
    }

    /// Build a `SnpLookup` the way `get_snp_position` does: `dsnp1` keyed
    /// `{reg1_pos}_{index}`, `dsnp2` keyed `{reg2_pos}_{index}`, both holding
    /// the same `{reg1_alleles}_{reg2_alleles}` value, and `dindex` mapping
    /// either key to the shared index.
    fn hap_lookup(sites: &[(i64, i64, &str)]) -> SnpLookup {
        let mut dsnp1 = indexmap::IndexMap::new();
        let mut dsnp2 = indexmap::IndexMap::new();
        let mut dindex: HashMap<String, usize> = HashMap::new();
        for (idx, (reg1_pos, reg2_pos, alleles)) in sites.iter().enumerate() {
            let k1 = format!("{reg1_pos}_{idx}");
            let k2 = format!("{reg2_pos}_{idx}");
            dsnp1.insert(k1.clone(), (*alleles).to_string());
            dsnp2.insert(k2.clone(), (*alleles).to_string());
            dindex.insert(k1, idx);
            dindex.insert(k2, idx);
        }
        SnpLookup {
            dsnp1,
            dsnp2,
            nchr: "chr22".to_string(),
            dindex,
        }
    }

    /// Scratch directory unique to the calling test, removed on drop.
    ///
    /// The previous `fn scratch(tag) -> PathBuf` was unlinked by the test body
    /// itself, which covered only the path where nothing panicked first: an
    /// `unwrap` in `write_hap_bam` or in `IndexedReader::from_path` leaked a
    /// BAM + BAI under `env::temp_dir()`. `Drop` runs during unwind too, so
    /// this covers every exit path, assertion failures included. Each test
    /// still drops its `bam::IndexedReader` explicitly before asserting, so
    /// the guard never unlinks a file the test holds open.
    struct Scratch(std::path::PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Self {
            let dir = std::env::temp_dir()
                .join(format!("cyrius_hap_{tag}_{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            Self(dir)
        }

        fn join(&self, name: &str) -> std::path::PathBuf {
            self.0.join(name)
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn a_read_whose_last_base_completes_the_site_keeps_its_haplotype() {
        // Review finding #7. `get_bases_per_read` skipped a read when
        // `end_pos >= seq.len()`, so a read whose extracted bases run to its
        // LAST stored base was dropped from the haplotypes -- even though
        // `seq[qpos..end_pos]` is a perfectly valid slice there, and even
        // though `snp_count`'s bound (`qpos + max_allele_len > seq.len()`)
        // still counts that same read in the D6/D7 depth ratios. The two
        // copies of this pileup loop had drifted apart on exactly this test.
        //
        // The fixture mirrors `call_var42127803hap`'s real panel: the two 1 bp
        // sites g.42127803 and g.42127941, 138 bp apart, read over 150 bp
        // reads. 138 bp apart in a 150 bp read leaves at most 11 bases of
        // slack past the second site, so a read covering both sites is always
        // close to its end there -- `second_site_at_read_end` starts so that
        // site 2 is its final stored base, i.e. `qpos + min_len == seq.len()`
        // and the bound has ZERO slack. It is the read that pins the `>`:
        // flipping it back to `>=` makes its haplotype "1x" instead of "12",
        // and since `extract_hap` discards any haplotype containing 'x' the
        // read then vanishes from the '12'/'21' counts the caller gates on.
        let dir = Scratch::new("end_bound");
        let bam_path = dir.join("hap_end_bound.bam");
        // Site 1 at 1-based 1001 (0-based 1000), site 2 at 1-based 1139
        // (0-based 1138) -- the same 138 bp spacing as 42127803/42127941.
        const READ_LEN: usize = 150;
        write_hap_bam(
            &bam_path,
            READ_LEN,
            &[
                // qpos 11 at site 1, qpos 149 at site 2: site 2 is this read's
                // final stored base.  1138 - 149 = 989.
                ("second_site_at_read_end", 989, &[(11, "C"), (149, "A")]),
                // qpos 0 / qpos 138: 11 bases of slack past site 2.
                ("roomy", 1000, &[(0, "T"), (138, "G")]),
            ],
        );

        let base_db = hap_lookup(&[(1001, 2001, "C_T"), (1139, 2139, "G_A")]);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1]).unwrap();
        let counts = extract_hap(&haps, &[0, 1]);

        // Close the BAM before asserting; `dir` unlinks itself on drop,
        // including when an assertion below panics.
        drop(reader);

        assert_eq!(
            haps.get("second_site_at_read_end").map(String::as_str),
            Some("12"),
            "site 2 is this read's last stored base, so `end_pos == seq.len()` \
             and the slice is still valid -- under the old `>=` bound the site \
             was dropped and the haplotype became \"1x\"; got {:?}",
            haps.get("second_site_at_read_end")
        );
        assert_eq!(
            haps.get("roomy").map(String::as_str),
            Some("21"),
            "control read with 11 bases of slack past site 2; got {:?}",
            haps.get("roomy")
        );
        // `call_var42127803hap` requires BOTH '12' and '21' to be seen.
        assert_eq!(
            counts.get("12").map(|v| v.iter().sum::<u32>()),
            Some(1),
            "the read ending exactly at site 2 must reach the '12' count; got {counts:?}"
        );
        assert_eq!(
            counts.get("21").map(|v| v.iter().sum::<u32>()),
            Some(1),
            "the roomy read must reach the '21' count; got {counts:?}"
        );
    }

    #[test]
    fn characterization_per_allele_length_fixes_commas_but_not_unequal_lengths() {
        // Review finding #7, second half. `min_len` used to come from the RAW
        // `reg1_allele`/`reg2_allele` strings, which can hold several
        // comma-separated alleles ('A,G'). The commas and the extra alleles
        // inflate the string length, so the single slice extracted here is
        // longer than the shortest allele at the site. Minimising over the
        // comma-SPLIT alleles of both regions is what makes the extraction
        // length agree with the matching.
        //
        // NOT "every read then scores 'x'" -- that overstatement was in the
        // production comment for a while and this very fixture refutes it. An
        // over-long slice misses every allele only when its length matches no
        // allele length at the site; when it happens to match a DIFFERENT
        // allele's length the read is scored for the WRONG region instead.
        // Both outcomes occur here, measured under the raw bound of 2:
        //   g_then_t -> "x"   (slice "GT" matches nothing)
        //   a_then_c -> "2"   (slice "AC" matches reg2's 2 bp allele)
        //
        // HONEST SCOPE: `data/CYP2D6_haplotype_38.txt` contains ZERO commas,
        // so this is a trap-removal test, not a bug reproduction. No shipped
        // input reaches it, and the 62-BAM before/after measurement cannot
        // move because of it.
        //
        // The fixture also shows the finding's own characterisation of when
        // the raw `min_len` goes wrong ("both alleles contain commas, or a
        // comma string is the shorter one") is too narrow: here only reg1 has
        // a comma and its raw string is the LONGER of the two ("A,G" is 3
        // bytes, "AC" is 2), yet the raw minimum is 2 while the true minimum
        // allele length is 1.
        let dir = Scratch::new("comma");
        let bam_path = dir.join("hap_comma.bam");
        const READ_LEN: usize = 50;
        write_hap_bam(
            &bam_path,
            READ_LEN,
            &[
                // Carries 'A' at the site, then 'C'.
                ("a_then_c", 980, &[(20, "AC")]),
                // Carries 'G' at the site, then 'T'.
                ("g_then_t", 990, &[(10, "GT")]),
            ],
        );

        // reg1 offers two 1 bp alleles, reg2 one 2 bp allele.
        //   raw  min_len = min("A,G".len(), "AC".len())        = 2  (wrong)
        //   split min_len = min(1 /*A*/, 1 /*G*/, 2 /*AC*/)    = 1  (right)
        let base_db = hap_lookup(&[(1001, 2001, "A,G_AC")]);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0]).unwrap();

        // Close the BAM before asserting; `dir` unlinks itself on drop,
        // including when an assertion below panics.
        drop(reader);

        assert_eq!(
            haps.get("g_then_t").map(String::as_str),
            Some("1"),
            "this read carries reg1's 1 bp 'G' allele. With the raw `min_len` \
             of 2 the extracted slice was \"GT\", which matches none of \
             reg1's 'A'/'G' or reg2's 'AC', so the read scored 'x' and \
             `extract_hap` threw it away; got {:?}",
            haps.get("g_then_t")
        );
        // RESIDUAL LIMITATION, pinned deliberately and NOT fixed here (out of
        // scope for this change): `get_bases_per_read` still extracts ONE
        // fixed-length slice per (read, site), so at a site whose alleles have
        // different lengths it can only ever test the shortest. This read
        // carries reg2's 2 bp "AC"; with the per-allele minimum of 1 the slice
        // is "A", which matches reg1's 1 bp 'A' instead, so the read now
        // scores '1' where the raw `min_len` of 2 happened to score it '2'.
        // Neither is correct -- a single fixed-length extraction cannot
        // discriminate alleles of unequal length at one site, which is what
        // `snp_count` does by matching each allele with its own length. This
        // is unreachable with today's data (see HONEST SCOPE above): the one
        // unequal-length row in the haplotype DB, line 10 / g.42130653
        // GCAGG vs GCAAGG, is consumed only by `call_var42130655ins_a` through
        // `get_supporting_reads_single_region` and never reaches
        // `get_bases_per_read`.
        assert_eq!(
            haps.get("a_then_c").map(String::as_str),
            Some("1"),
            "CHARACTERIZATION, KNOWN-WRONG -- this assertion is a \
             change-detector, NOT a specification. The read carries reg2's 2 bp \
             \"AC\" and is labelled region 1, which is wrong; the expected \
             value is pinned only so the residual unequal-length limitation \
             cannot change unnoticed. When per-allele matching lands (the \
             correct fix), DELETE this assertion rather than chasing it; got \
             {:?}",
            haps.get("a_then_c")
        );
    }
    #[test]
    fn a_trailing_comma_must_not_fabricate_a_match_on_every_read() {
        // Review of the `min_len` change: with the comma SPLIT, `min_len == 0`
        // became reachable, and it fails OPEN. One trailing comma in an allele
        // field -- `"A,"` splits to `["A", ""]` -- makes the minimum 0, so the
        // extracted slice is the EMPTY string. `get_hap_counts` splits
        // `base1`/`base2` on ',' as well, so it compares that empty slice
        // against an empty component and `"" == ""` MATCHES: without the
        // `min_len == 0` guard EVERY covering read at the site is scored '1'
        // and counted by `extract_hap`. The failure mode is a fabricated
        // haplotype on every read, not a discarded one.
        //
        // This exposure is strictly NEW, not inherited from the raw
        // `min(|reg1|, |reg2|)`: that needed an empty RAW field, which
        // `get_snp_position` cannot produce (its fields come from
        // `split_whitespace`). The split needs only an empty component.
        //
        // Measured without the guard (`min_len == 0 { continue }` removed):
        //   haps   = {"read_ac": "1", "read_gt": "1"}
        //   counts = {"1": [1, 1]}
        // i.e. 2 of 2 covering reads fabricated.
        let dir = Scratch::new("zero_len");
        let bam_path = dir.join("hap_zero_len.bam");
        const READ_LEN: usize = 50;
        write_hap_bam(
            &bam_path,
            READ_LEN,
            &[
                ("read_ac", 980, &[(20, "AC")]),
                ("read_gt", 990, &[(10, "GT")]),
            ],
        );

        // reg1 = "A," -> ["A", ""], reg2 = "AC" -> ["AC"]; the minimum over
        // those three components is 0.
        let base_db = hap_lookup(&[(1001, 2001, "A,_AC")]);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0]).unwrap();
        let counts = extract_hap(&haps, &[0]);

        // Close the BAM before asserting; `dir` unlinks itself on drop,
        // including when an assertion below panics.
        drop(reader);

        assert!(
            haps.is_empty(),
            "a zero-length extraction must be skipped, not compared: the empty \
             slice matches an empty comma component in `get_hap_counts` and \
             fabricates a '1' for every covering read. Without the \
             `min_len == 0` guard this is {{\"read_ac\": \"1\", \
             \"read_gt\": \"1\"}}; got {haps:?}"
        );
        assert!(
            counts.is_empty(),
            "no haplotype may reach `extract_hap` from a zero-length slice; \
             without the guard this is {{\"1\": [1, 1]}}, i.e. 2 of 2 covering \
             reads counted on nothing; got {counts:?}"
        );
    }

    #[test]
    fn a_multi_base_site_at_the_read_end_keeps_its_haplotype() {
        // Same `>=` -> `>` bound as
        // `a_read_whose_last_base_completes_the_site_keeps_its_haplotype`, but
        // at MULTI-BASE sites, which that 1 bp fixture does not exercise at
        // all: the slice there is `seq[qpos..qpos + 1]`, a single base, so it
        // says nothing about a 10 or 12 base slice landing flush with the end
        // of the stored sequence.
        //
        // Exactly one reachable panel has multi-base rows -- `g.42126938C>T`,
        // `data/CYP2D6_haplotype_38.txt` lines 2 and 3, 10 bp and 12 bp -- and
        // it is well covered on real data. Replaying the three reachable panels
        // over the 62 benchmark BAMs, 99 of the 261 newly admitted observations
        // (38%) sit at those two 10/12 bp sites (plus their D7 mates), spread
        // over 47 of the 62 samples; HG006 contributes 6, HG02373 4 and
        // NA19777 1.
        //
        // Note what is NOT true: multi-base sites are not MORE exposed per
        // site. Both bounds reject `qpos + min_len > seq.len()`, so they differ
        // on the single value `qpos == seq.len() - min_len` whatever `min_len`
        // is, and the measured rescue count per (site, position) is about the
        // same at 1 bp and at 10/12 bp sites (20-34 reads each). What is
        // different is the slice itself, and that is what this test covers.
        //
        // Geometry copied from that panel: a 12 bp site and a 10 bp site 47 bp
        // apart (as g.42126611 and g.42126658 are), read over 150 bp reads.
        // `site_at_read_end` is placed so the 10 bp slice ends at exactly
        // `seq.len()` -- zero slack. Under `>=` its second site is dropped and
        // the haplotype is "1x", which `extract_hap` then discards.
        let dir = Scratch::new("multibase_end");
        let bam_path = dir.join("hap_multibase_end.bam");
        const READ_LEN: usize = 150;
        // Site A (index 0) at 1-based 1001 = 0-based 1000, 12 bp alleles.
        // Site B (index 1) at 1-based 1048 = 0-based 1047, 10 bp alleles.
        const A_REG1: &str = "GTCACCAGGAAA"; // 12 bp, as line 3 of the DB
        const A_REG2: &str = "CTCACCAGGAAA";
        const B_REG1: &str = "AGTGGGCACC"; // 10 bp, as line 2 of the DB
        const B_REG2: &str = "GGCGGCCACG";
        write_hap_bam(
            &bam_path,
            READ_LEN,
            &[
                // 0-based start 907: site A at qpos 93, site B at qpos 140, so
                // 140 + 10 == 150 == seq.len() and the bound has ZERO slack.
                (
                    "site_at_read_end",
                    907,
                    &[(93, A_REG1), (140, B_REG2)],
                ),
                // 0-based start 1000: site A at qpos 0, site B at qpos 47, so
                // 93 bases of slack past the 10 bp site.
                ("roomy", 1000, &[(0, A_REG2), (47, B_REG1)]),
            ],
        );

        let base_db = hap_lookup(&[
            (1001, 2001, "GTCACCAGGAAA_CTCACCAGGAAA"),
            (1048, 2048, "AGTGGGCACC_GGCGGCCACG"),
        ]);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0, 1]).unwrap();
        let counts = extract_hap(&haps, &[0, 1]);

        // Close the BAM before asserting; `dir` unlinks itself on drop,
        // including when an assertion below panics.
        drop(reader);

        assert_eq!(
            haps.get("site_at_read_end").map(String::as_str),
            Some("12"),
            "the 10 bp site is this read's final 10 stored bases, so \
             `end_pos == seq.len()` and `seq[qpos..end_pos]` is a valid slice; \
             under the old `>=` the site was dropped and the haplotype became \
             \"1x\"; got {:?}",
            haps.get("site_at_read_end")
        );
        assert_eq!(
            haps.get("roomy").map(String::as_str),
            Some("21"),
            "control read with 93 bases of slack past the 10 bp site; got {:?}",
            haps.get("roomy")
        );
        assert_eq!(
            counts.get("12").map(|v| v.iter().sum::<u32>()),
            Some(1),
            "the read ending exactly at the 10 bp site must reach the '12' \
             count; got {counts:?}"
        );
        assert_eq!(
            counts.get("21").map(|v| v.iter().sum::<u32>()),
            Some(1),
            "the roomy read must reach the '21' count; got {counts:?}"
        );
    }

    #[test]
    fn both_regions_comma_separated_with_equal_length_alleles_are_still_scored() {
        // The case that makes the comma SPLIT strictly better than the raw
        // `min(|reg1_allele|, |reg2_allele|)`, and the justification for
        // keeping that half of the change.
        //
        // When both regions are comma'd and every allele has the same length,
        // the raw minimum matches NO allele length at all: here 'A,G' and
        // 'C,T' are both 3 bytes, so `raw == 3` while every allele is 1 bp. The
        // 3 bp slice can never equal a 1 bp allele, so every read at the site
        // scores 'x', `extract_hap` returns EMPTY, the whole panel reads as "no
        // reads", and every gate built on it defaults false -- a silent, total
        // loss of the site rather than a marginal miscount.
        //
        // Measured by running this function under the raw bound:
        //   haps   = {"carries_a": "x", "carries_t": "x"}
        //   counts = {}
        // and under the split bound (asserted below):
        //   haps   = {"carries_a": "1", "carries_t": "2"}
        //   counts = {"1": [1], "2": [1]}
        //
        // HONEST SCOPE: ZERO of the 285 rows across the four shipped position
        // DBs have both regions comma'd, and the haplotype DB this function
        // reads has no commas at all, so this is trap removal, not a bug
        // reproduction. No benchmark number can move because of it.
        let dir = Scratch::new("both_comma");
        let bam_path = dir.join("hap_both_comma.bam");
        const READ_LEN: usize = 50;
        write_hap_bam(
            &bam_path,
            READ_LEN,
            &[
                // Only the site base is set; the rest of the read is 'N', so a
                // 3 bp slice reads "ANN" / "TNN" and matches nothing.
                ("carries_a", 980, &[(20, "A")]),
                ("carries_t", 990, &[(10, "T")]),
            ],
        );

        //   raw   min_len = min("A,G".len(), "C,T".len()) = 3  (matches nothing)
        //   split min_len = min(1, 1, 1, 1)               = 1  (right)
        let base_db = hap_lookup(&[(1001, 2001, "A,G_C,T")]);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let haps = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0]).unwrap();
        let counts = extract_hap(&haps, &[0]);

        // Close the BAM before asserting; `dir` unlinks itself on drop,
        // including when an assertion below panics.
        drop(reader);

        assert_eq!(
            haps.get("carries_a").map(String::as_str),
            Some("1"),
            "this read carries reg1's 1 bp 'A'; under the raw 3 bp bound the \
             slice was \"ANN\" and the read scored 'x'; got {:?}",
            haps.get("carries_a")
        );
        assert_eq!(
            haps.get("carries_t").map(String::as_str),
            Some("2"),
            "this read carries reg2's 1 bp 'T'; under the raw 3 bp bound the \
             slice was \"TNN\" and the read scored 'x'; got {:?}",
            haps.get("carries_t")
        );
        assert_eq!(
            counts.get("1").map(|v| v.iter().sum::<u32>()),
            Some(1),
            "under the raw bound `extract_hap` was EMPTY here -- the entire \
             site read as 'no reads'; got {counts:?}"
        );
        assert_eq!(
            counts.get("2").map(|v| v.iter().sum::<u32>()),
            Some(1),
            "under the raw bound `extract_hap` was EMPTY here -- the entire \
             site read as 'no reads'; got {counts:?}"
        );
    }

    // ---- E1: a failed BAM fetch is a no-call, not zero reads ----
    //
    // Escalation E1 of review batch 4, decided by the human 2026-10-07: when
    // htslib cannot position the reader (broken or stale index, CRAM reference
    // that cannot be resolved, invalid interval) the reads were never looked
    // at. Returning zero counts makes that look like evidence, and the sample
    // gets a confident call. It must surface as an Err, which main.rs turns
    // into a no-call. An inverted interval is the input used here because it
    // is the one that reliably drives `fetch` into Err through the public API
    // (measured: Region(tid, 100, 50) errs; an out-of-range tid or a negative
    // start does not, they just return no reads).

    #[test]
    fn get_haplotypes_is_an_error_when_the_fetch_fails() {
        // As in snp_count: a site at -5 makes the pileup interval inverted.
        let dir = Scratch::new("fetch_err");
        let bam_path = dir.join("fetch_err.bam");
        write_hap_bam(&bam_path, 150, &[("r1", 1000, &[(0, "C")])]);
        let base_db = hap_lookup(&[(-5, 2001, "C_T")]);
        let mut reader = bam::IndexedReader::from_path(&bam_path).unwrap();

        let r = get_haplotypes_from_bam_single_region(&mut reader, &base_db, &[0]);
        drop(reader);

        assert!(r.is_err(), "a failed fetch must be an error, not an empty haplotype map");
    }
}

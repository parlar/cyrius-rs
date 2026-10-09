//! Phase audit: the read counts behind `--phase-disambiguate`, from the
//! shipped counting code, for a BAM and a list of variant pairs or one site.
//!
//! Usage:
//!   phase_audit pairs <bam> <pairs.tsv>
//!     `pairs.tsv` has a header; the columns named `va` and `vb` hold variant
//!     names as in the target variant file (column 7). Other columns are
//!     ignored. Prints one TSV row per pair:
//!     va vb distance aa ar ra rr_all rr_anchored cis trans
//!     (`NA` in every count column when a variant has no lookup entry).
//!   phase_audit site <bam> <variant_name>
//!     Prints `variant alt ref none`: reads overlapping the site called ALT,
//!     REF, or neither.
//!
//! Pairs are counted on the single-fetch path `disambiguate` uses by default
//! (no `--phase-readpair`). The contig is `chr22` if the header has it, else
//! `22`.

use cyrius_rs::caller::phase_disambiguate::{count_pair, count_site, VariantLookup};
use cyrius_rs::data;
use cyrius_rs::depth_calling::utilities::open_alignment_file_with_index;
use rust_htslib::bam::{self, Read};
use std::process::exit;

const USAGE: &str = "usage: phase_audit pairs <bam> <pairs.tsv>\n       phase_audit site <bam> <variant_name>";

fn die(message: &str) -> ! {
    eprintln!("phase_audit: {message}");
    exit(2);
}

fn open(bam_path: &str) -> (bam::IndexedReader, &'static str) {
    let reader = open_alignment_file_with_index(bam_path, None, None)
        .unwrap_or_else(|e| die(&format!("cannot open {bam_path}: {e}")));
    let nchr = if reader.header().tid(b"chr22").is_some() { "chr22" } else { "22" };
    (reader, nchr)
}

fn audit_pairs(bam_path: &str, pairs_path: &str, lookup: &VariantLookup) {
    let text = std::fs::read_to_string(pairs_path)
        .unwrap_or_else(|e| die(&format!("cannot read {pairs_path}: {e}")));
    let mut lines = text.lines();
    let header: Vec<&str> = lines
        .next()
        .unwrap_or_else(|| die(&format!("{pairs_path} is empty")))
        .split('\t')
        .collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|&c| c == name)
            .unwrap_or_else(|| die(&format!("{pairs_path} has no '{name}' column")))
    };
    let (ia, ib) = (column("va"), column("vb"));

    let (mut reader, nchr) = open(bam_path);
    println!("va\tvb\tdistance\taa\tar\tra\trr_all\trr_anchored\tcis\ttrans");
    for (n, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (va, vb) = match (fields.get(ia), fields.get(ib)) {
            (Some(a), Some(b)) => (*a, *b),
            _ => die(&format!("{pairs_path} line {}: too few columns", n + 2)),
        };
        match count_pair(&mut reader, nchr, va, vb, lookup) {
            Some(e) => {
                let c = e.counts;
                println!(
                    "{va}\t{vb}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    e.distance_bp(),
                    c.aa,
                    c.ar,
                    c.ra,
                    c.rr_all,
                    c.rr_anchored,
                    c.cis(),
                    c.trans()
                );
            }
            None => println!("{va}\t{vb}{}", "\tNA".repeat(8)),
        }
    }
}

fn audit_site(bam_path: &str, variant: &str, lookup: &VariantLookup) {
    let (mut reader, nchr) = open(bam_path);
    println!("variant\talt\tref\tnone");
    match count_site(&mut reader, nchr, variant, lookup) {
        Some((alt, ref_count, none)) => println!("{variant}\t{alt}\t{ref_count}\t{none}"),
        None => println!("{variant}{}", "\tNA".repeat(3)),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let lookup = VariantLookup::from_target_file(data::TARGET_VARIANT);
    match args.get(1).map(String::as_str) {
        Some("pairs") if args.len() == 4 => audit_pairs(&args[2], &args[3], &lookup),
        Some("site") if args.len() == 4 => audit_site(&args[2], &args[3], &lookup),
        _ => die(USAGE),
    }
}

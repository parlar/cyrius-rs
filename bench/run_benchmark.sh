#!/usr/bin/env bash
# Run cyrius-rs on every GeT-RM sample in the truth file and score the calls.
#
# Usage:
#   bash bench/run_benchmark.sh [FLAGS...]
#
# Environment:
#   BENCH_BAM_DIR  directory with {sample}_pgx_regions.bam (default: <repo>/benchmark)
#   BENCH_JOBS     parallel jobs (default: nproc)
#   BENCH_KEEP=1   keep the per-sample JSON/logs and print where they are
#
# Each sample is scored as one of:
#   EXACT      a single genotype that matches the truth
#   AMBIGUOUS  several ';' genotypes, one of which matches the truth
#   WRONG      no reported genotype matches the truth
#   NO_CALL    the tool ran but reported no genotype
#   ERROR      the tool exited non-zero or wrote no readable JSON
#   NO_BAM     no BAM for this truth sample in BENCH_BAM_DIR
# EXACT is the headline number; AMBIGUOUS is reported separately.
#
# Before comparing, both sides are normalized: haplotypes are sorted,
# "*A(*B)" parent notation becomes "*B" (truth for NA19143 is written as
# *2(*45)), and "*AxN" becomes N copies joined by "+".
# Exit status is 1 if any sample is ERROR, else 0.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARY="$PROJECT_ROOT/cyrius-rs/target/release/cyrius-rs"
TRUTH_FILE="$SCRIPT_DIR/getrm_cyp2d6_truth.tsv"
BAM_DIR="${BENCH_BAM_DIR:-$PROJECT_ROOT/benchmark}"
JOBS="${BENCH_JOBS:-$(nproc)}"
RESULTS_DIR=$(mktemp -d "${TMPDIR:-/tmp}/bench_XXXXXX")
if [ "${BENCH_KEEP:-0}" != "1" ]; then
    trap 'rm -rf "$RESULTS_DIR"' EXIT
fi

if [ ! -f "$BINARY" ]; then
    echo "Binary not found. Building..."
    (cd "$PROJECT_ROOT/cyrius-rs" && cargo build --release 2>&1 | tail -1)
fi

mapfile -t SAMPLES < <(grep -v '^#' "$TRUTH_FILE" | cut -f1 | grep -v '^$')
echo "Truth samples: ${#SAMPLES[@]}   BAM dir: $BAM_DIR   jobs: $JOBS"
if [ $# -gt 0 ]; then
    echo "Flags: $*"
fi
echo ""

printf '%s\n' "${SAMPLES[@]}" | xargs -P "$JOBS" -I{} bash -c '
    sample="$1"; shift
    binary="$1"; shift
    bam_dir="$1"; shift
    outdir="$1"; shift
    bam="$bam_dir/${sample}_pgx_regions.bam"
    if [ ! -f "$bam" ]; then
        echo "NO_BAM" > "$outdir/$sample.status"
        exit 0
    fi
    "$binary" -i "$bam" -o "$outdir" --id "$sample" "$@" 2> "$outdir/$sample.log"
    echo "$?" > "$outdir/$sample.status"
' _ {} "$BINARY" "$BAM_DIR" "$RESULTS_DIR" "$@" || true

status=0
python3 - "$RESULTS_DIR" "$TRUTH_FILE" <<'EOF' || status=$?
import json, os, re, sys

results_dir, truth_file = sys.argv[1:]

def normalize(diplotype):
    s = re.sub(r"\*[0-9.]+\((\*[0-9.]+)\)", r"\1", diplotype.strip())
    s = re.sub(r"(\*[0-9.]+)x(\d+)", lambda m: "+".join([m.group(1)] * int(m.group(2))), s)
    return "/".join(sorted(s.split("/")))

truth = {}
for line in open(truth_file):
    if line.startswith("#") or not line.strip():
        continue
    fields = line.rstrip("\n").split("\t")
    truth[fields[0]] = fields[1]

rows = []
for sample, expected in truth.items():
    status_path = os.path.join(results_dir, sample + ".status")
    status = open(status_path).read().strip() if os.path.exists(status_path) else "1"
    called = None
    if status == "NO_BAM":
        category = "NO_BAM"
    else:
        try:
            data = json.load(open(os.path.join(results_dir, sample + ".json")))
            called = data[sample].get("Genotype")
            readable = True
        except (OSError, ValueError, KeyError):
            readable = False
        if status != "0" or not readable:
            category = "ERROR"
        elif not called:
            category = "NO_CALL"
        else:
            options = [normalize(o) for o in called.split(";")]
            if normalize(expected) not in options:
                category = "WRONG"
            elif len(options) == 1:
                category = "EXACT"
            else:
                category = "AMBIGUOUS"
    rows.append((sample, called, expected, category))

print(f"{'Sample':<12} {'Called':<34} {'Truth':<26} Status")
print(f"{'------':<12} {'------':<34} {'-----':<26} ------")
for sample, called, expected, category in rows:
    if category != "EXACT":
        print(f"{sample:<12} {str(called):<34} {expected:<26} {category}")

counts = {c: sum(1 for r in rows if r[3] == c) for c in
          ["EXACT", "AMBIGUOUS", "WRONG", "NO_CALL", "ERROR", "NO_BAM"]}
scored = len(rows) - counts["NO_BAM"]
print()
print("=== Summary ===")
print(f"Exact:      {counts['EXACT']} / {scored}")
print(f"Ambiguous:  {counts['AMBIGUOUS']}  (truth is one of several reported genotypes)")
print(f"Wrong:      {counts['WRONG']}")
print(f"No call:    {counts['NO_CALL']}")
print(f"Error:      {counts['ERROR']}")
if counts["NO_BAM"]:
    print(f"No BAM:     {counts['NO_BAM']}")
sys.exit(1 if counts["ERROR"] else 0)
EOF

if [ "${BENCH_KEEP:-0}" = "1" ]; then
    echo "Results kept in $RESULTS_DIR"
fi
exit $status

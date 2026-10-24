from __future__ import annotations

import argparse
import csv
import gzip
import sys
import time
from pathlib import Path

OUT_MIN_DP = 10
OUT_MIN_GQ = 20
MID_GQ_MAX = 50
OUT_MIN_CADD = 20.0
WIN_MIN_CADD = 10.0
AB_LOW = 0.2
AB_HIGH = 0.8
OUTCOMES = ("kept", "kept_positive_control", "no_cadd_list", "removed_quality", "removed_cadd")


def header_lines(controls: str) -> list[str]:
    return [
        '##INFO=<ID=F2_FLAGS,Number=.,Type=String,Description="Filter 2 flags on kept records: NOT_PASS, '
        f'LOW_DP (DP<={OUT_MIN_DP}), LOW_GQ (GQ<={OUT_MIN_GQ}), MID_GQ (GQ {OUT_MIN_GQ + 1}-{MID_GQ_MAX}), AB_SKEW (het allele balance outside {AB_LOW}-{AB_HIGH}), '
        'NO_CADD, CADD_BELOW_RULE, POSITIVE_CONTROL">\n',
        '##INFO=<ID=F2_PC,Number=1,Type=String,Description="Positive control ID (docs/positive_controls_grouped.tsv) '
        'matched by exact CHROM/POS/REF/ALT">\n',
        f'##Filter2=<Version="2.0",Outside="DP>{OUT_MIN_DP},GQ>{OUT_MIN_GQ},CADD_PHRED>={OUT_MIN_CADD:g}",'
        f'MidGQ="flag {OUT_MIN_GQ + 1}-{MID_GQ_MAX}",MvaWindow="CADD_PHRED>={WIN_MIN_CADD:g},DP/GQ flag only",PASS="flag only",NoCADD="separate file",'
        f'PositiveControls="{controls} always kept",Spec="docs/14_FILTER2_QUALITY_CADD.md">\n',
    ]


def to_num(v: str | None) -> float | None:
    return float(v) if v not in (None, "", ".") else None


def load_controls(path: Path) -> dict[tuple[str, str, str, str], str]:
    with open(path, encoding="utf-8", newline="") as fh:
        return {(r["chrom"], r["pos"], r["ref"], r["alt"]): r["id"] for r in csv.DictReader(fh, delimiter="\t")}


def quality_flags(cols: list[str]) -> tuple[list[str], bool]:
    sample = dict(zip(cols[8].split(":"), cols[9].split(":")))
    gt = sample.get("GT", ".").replace("|", "/")
    dp, gq = to_num(sample.get("DP")), to_num(sample.get("GQ"))
    flags = []
    if cols[6] != "PASS":
        flags.append("NOT_PASS")
    if dp is None or dp <= OUT_MIN_DP:
        flags.append("LOW_DP")
    if gq is None or gq <= OUT_MIN_GQ:
        flags.append("LOW_GQ")
    elif gq <= MID_GQ_MAX:
        flags.append("MID_GQ")
    if gt in ("0/1", "1/0"):
        ad = [to_num(x) for x in sample.get("AD", ".").split(",")]
        if len(ad) == 2 and None not in ad and ad[0] + ad[1] > 0:
            ab = ad[1] / (ad[0] + ad[1])
            if ab < AB_LOW or ab > AB_HIGH:
                flags.append("AB_SKEW")
    return flags, "LOW_DP" in flags or "LOW_GQ" in flags


def classify(cols: list[str], controls: dict) -> tuple[str, list[str], str | None]:
    info = cols[7].split(";")
    f1 = next((x[10:] for x in info if x.startswith("F1_REASON=")), None)
    if f1 is None:
        sys.exit(f"ERROR: record without F1_REASON at {cols[0]}:{cols[1]}")
    cadd = to_num(next((x[11:] for x in info if x.startswith("CADD_PHRED=")), None))
    in_window = "MVA_WINDOW" in f1.split(",")
    flags, low_dp_gq = quality_flags(cols)
    if cadd is None:
        flags.append("NO_CADD")
    elif cadd < (WIN_MIN_CADD if in_window else OUT_MIN_CADD):
        flags.append("CADD_BELOW_RULE")
    pc = controls.get((cols[0], cols[1], cols[3], cols[4]))
    if pc:
        flags.append("POSITIVE_CONTROL")
        return "kept_positive_control", flags, pc
    if low_dp_gq and not in_window:
        return "removed_quality", flags, None
    if cadd is None:
        return "no_cadd_list", flags, None
    if "CADD_BELOW_RULE" in flags:
        return "removed_cadd", flags, None
    return "kept", flags, None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="CADD-annotated normalized Filter 1 VCF (.vcf or .vcf.gz)")
    ap.add_argument("output", help="main output VCF (uncompressed)")
    ap.add_argument("no_cadd_output", help="VCF for records that pass quality but have no CADD score")
    ap.add_argument("--controls", default=str(Path(__file__).resolve().parents[1] / "config" / "positive_controls_grouped.tsv"))
    args = ap.parse_args()

    t0 = time.time()
    controls = load_controls(Path(args.controls))
    opener = gzip.open if args.input.endswith(".gz") else open
    counts = {f"{scope}_{o}": 0 for scope in ("outside", "window") for o in OUTCOMES}
    n_in = 0
    saw_cadd_header = False

    with opener(args.input, "rt", encoding="utf-8", errors="surrogateescape", newline="") as fin, \
         open(args.output, "w", encoding="utf-8", errors="surrogateescape", newline="") as fmain, \
         open(args.no_cadd_output, "w", encoding="utf-8", errors="surrogateescape", newline="") as fnone:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("##INFO=<ID=CADD_PHRED,"):
                    saw_cadd_header = True
                if line.startswith("#CHROM"):
                    if not saw_cadd_header:
                        sys.exit("ERROR: no ##INFO=<ID=CADD_PHRED> header; run scripts/cadd_annotate.sh first")
                    for f in (fmain, fnone):
                        f.writelines(header_lines(Path(args.controls).name))
                for f in (fmain, fnone):
                    f.write(line)
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            outcome, flags, pc = classify(cols, controls)
            scope = "window" if "MVA_WINDOW" in line else "outside"
            counts[f"{scope}_{outcome}"] += 1
            if outcome.startswith("removed"):
                continue
            if flags:
                cols[7] += ";F2_FLAGS=" + ",".join(flags)
            if pc:
                cols[7] += ";F2_PC=" + pc
            (fnone if outcome == "no_cadd_list" else fmain).write("\t".join(cols) + "\n")

    main_n = sum(v for k, v in counts.items() if k.endswith(("_kept", "_kept_positive_control")))
    none_n = sum(v for k, v in counts.items() if k.endswith("_no_cadd_list"))
    summary = [f"records_in\t{n_in}", f"main_output\t{main_n}", f"no_cadd_output\t{none_n}",
               f"removed\t{n_in - main_n - none_n}", f"positive_controls_loaded\t{len(controls)}"]
    summary += [f"{k}\t{v}" for k, v in counts.items()]
    print("\n".join(summary), file=sys.stderr)
    print(f"elapsed_seconds\t{time.time() - t0:.0f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

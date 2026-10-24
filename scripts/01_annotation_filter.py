
from __future__ import annotations

import argparse
import gzip
import re
import sys
import time
from pathlib import Path

SPLICE_TERMS = frozenset({
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
})
PLP_WORDS = frozenset({"pathogenic", "likely_pathogenic"})
CONFLICTING_WORD = "conflicting_classifications_of_pathogenicity"
REASONS = ("HIGH_MODERATE", "SPLICE", "MVA_WINDOW", "CLINVAR_PLP", "CLINVAR_CONFLICTING")
CLINSIG_SPLIT = re.compile(r"[&/]")


def header_lines(window: int) -> list[str]:
    return [
        '##INFO=<ID=F1_REASON,Number=.,Type=String,Description="Filter 1 (annotation) rules '
        'that kept this record: HIGH_MODERATE, SPLICE, MVA_WINDOW, CLINVAR_PLP, CLINVAR_CONFLICTING">\n',
        f'##Filter1=<Version="1.0",Window={window},Spec="docs/reproducibility.md",'
        'CSQ="unchanged",QualityRarity="not applied">\n',
    ]


def load_windows(path: Path, window: int) -> dict[str, list[tuple[int, int]]]:
    windows: dict[str, list[tuple[int, int]]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        gene, chrom, start, end = line.split("\t")[:4]
        windows.setdefault(chrom, []).append((int(start) - window, int(end) + window))
    return windows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="VEP-annotated VCF (.vcf or .vcf.gz)")
    ap.add_argument("output", help="output VCF (uncompressed)")
    ap.add_argument("--genes", default=str(Path(__file__).resolve().parents[1] / "config" / "mva_genes.tsv"))
    ap.add_argument("--window", type=int, default=100000)
    args = ap.parse_args()

    t0 = time.time()
    windows = load_windows(Path(args.genes), args.window)
    opener = gzip.open if args.input.endswith(".gz") else open
    counts = dict.fromkeys(REASONS, 0)
    n_in = n_out = 0
    idx = None

    with opener(args.input, "rt", encoding="utf-8", errors="surrogateescape", newline="") as fin, \
         open(args.output, "w", encoding="utf-8", errors="surrogateescape", newline="") as fout:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("##INFO=<ID=CSQ,"):
                    fmt = line.split("Format: ", 1)[1].split('">', 1)[0].split("|")
                    idx = {name: i for i, name in enumerate(fmt)}
                elif line.startswith("#CHROM"):
                    if idx is None:
                        sys.exit("ERROR: no ##INFO=<ID=CSQ> header found")
                    fout.writelines(header_lines(args.window))
                fout.write(line)
                continue

            n_in += 1
            tab1 = line.index("\t")
            tab2 = line.index("\t", tab1 + 1)
            chrom, pos = line[:tab1], int(line[tab1 + 1:tab2])
            in_window = any(lo <= pos <= hi for lo, hi in windows.get(chrom, ()))

            if not in_window and not ("HIGH" in line or "MODERATE" in line
                                      or "splice_" in line or "pathogenic" in line):
                continue

            cols = line.rstrip("\n").split("\t")
            info = cols[7].split(";")
            csq = next((f[4:] for f in info if f.startswith("CSQ=")), "")
            hm = sp = plp = conf = False
            if csq:
                i_imp, i_can = idx["IMPACT"], idx["CANONICAL"]
                i_csq, i_sig = idx["Consequence"], idx["CLIN_SIG"]
                for ann in csq.split(","):
                    a = ann.split("|")
                    if a[i_can] == "YES":
                        if a[i_imp] in ("HIGH", "MODERATE"):
                            hm = True
                        if not sp and not SPLICE_TERMS.isdisjoint(a[i_csq].split("&")):
                            sp = True
                    if a[i_sig]:
                        words = CLINSIG_SPLIT.split(a[i_sig])
                        if not plp and not PLP_WORDS.isdisjoint(words):
                            plp = True
                        if not conf and CONFLICTING_WORD in words:
                            conf = True

            reasons = [r for r, hit in zip(REASONS, (hm, sp, in_window, plp, conf)) if hit]
            if not reasons:
                continue
            for r in reasons:
                counts[r] += 1
            cols[7] += ";F1_REASON=" + ",".join(reasons)
            fout.write("\t".join(cols) + "\n")
            n_out += 1

    summary = [f"records_in\t{n_in}", f"records_kept\t{n_out}"]
    summary += [f"reason_{r}\t{counts[r]}" for r in REASONS]
    print("\n".join(summary), file=sys.stderr)
    print(f"elapsed_seconds\t{time.time() - t0:.0f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())


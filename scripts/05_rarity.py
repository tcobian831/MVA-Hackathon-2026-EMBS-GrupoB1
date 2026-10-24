from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from quality_helpers import load_controls, to_num

MAX_AF = 0.01
GRPMAX_GROUPS = ("AFR", "AMR", "EAS", "MID", "NFE", "SAS")
EXCLUDED_GROUPS = ("AMI", "ASJ", "FIN", "REMAINING")
BUB1B_SENSITIVITY_AF = 0.05
OUTCOMES = ("kept", "kept_no_freq", "kept_positive_control", "removed_common")


def header_lines(controls: str) -> list[str]:
    groups = ",".join(GRPMAX_GROUPS)
    return [
        f'##INFO=<ID=F3_GRPMAX,Number=1,Type=Float,Description="Highest gnomAD v4.1 genomes AF of this ALT among {groups} '
        '(the gnomAD v4.1 grpmax groups), from the CSQ entries whose Allele matches this ALT">\n',
        '##INFO=<ID=F3_GRPMAX_POP,Number=1,Type=String,Description="Group that gives F3_GRPMAX">\n',
        '##INFO=<ID=F3_FLAGS,Number=.,Type=String,Description="Filter 3 flags on kept records: NO_FREQ (no gnomAD genomes AF at all), '
        f'NO_GRPMAX (AF only in excluded groups), COMMON_IN_EXCLUDED_GROUP (AF>={MAX_AF:g} in {",".join(EXCLUDED_GROUPS)}), '
        f'GLOBAL_AF_COMMON (gnomADg_AF>={MAX_AF:g}), POSITIVE_CONTROL">\n',
        f'##Filter3=<Version="1.0",Rule="F3_GRPMAX<{MAX_AF:g}",Missing="kept",Groups="{groups}",Excluded="{",".join(EXCLUDED_GROUPS)}",'
        f'AlleleMatched="yes",PositiveControls="{controls} always kept",Spec="docs/15_FILTER3_RARITY.md">\n',
    ]


def vep_allele(ref: str, alts: list[str], idx: int) -> str:
    alleles = [ref] + alts
    if len({a[0] for a in alleles}) == 1:
        alleles = [a[1:] or "-" for a in alleles]
    return alleles[idx]


def own_allele(cols: list[str], info: dict[str, str]) -> str:
    orig = info.get("ORIGINAL")
    if orig:
        _, _, ref, alts, idx = orig.split("|")
        return vep_allele(ref, alts.split(","), int(idx))
    return vep_allele(cols[3], [cols[4]], 1)


def frequencies(cols: list[str], info: dict[str, str], names: list[str]) -> tuple[dict[str, float | None], set[str]]:
    allele = own_allele(cols, info)
    fields = ["gnomADg_AF"] + [f"gnomADg_{g}_AF" for g in GRPMAX_GROUPS + EXCLUDED_GROUPS]
    idx = {f: names.index(f) for f in fields}
    a_idx, s_idx = names.index("Allele"), names.index("SYMBOL")
    seen, symbols = set(), set()
    for entry in info["CSQ"].split(","):
        f = entry.split("|")
        if f[a_idx] != allele:
            continue
        seen.add(tuple(f[idx[k]] for k in fields))
        symbols.add(f[s_idx])
    where = f"{cols[0]}:{cols[1]} {cols[3]}>{cols[4]}"
    if not seen:
        sys.exit(f"ERROR: no CSQ entry with Allele={allele} at {where}")
    if len(seen) > 1:
        sys.exit(f"ERROR: CSQ entries for Allele={allele} disagree on gnomAD values at {where}")
    values = seen.pop()
    if any("&" in v for v in values):
        sys.exit(f"ERROR: multi-valued gnomAD field at {where}")
    return {k: to_num(v) for k, v in zip(fields, values)}, symbols


def classify(cols: list[str], names: list[str], controls: dict) -> tuple[str, list[str], float | None, str | None, dict]:
    info = dict(x.split("=", 1) if "=" in x else (x, "") for x in cols[7].split(";"))
    af, symbols = frequencies(cols, info, names)
    grp = {g: af[f"gnomADg_{g}_AF"] for g in GRPMAX_GROUPS if af[f"gnomADg_{g}_AF"] is not None}
    pop = max(grp, key=grp.get) if grp else None
    grpmax = grp[pop] if pop else None
    excl = [af[f"gnomADg_{g}_AF"] for g in EXCLUDED_GROUPS if af[f"gnomADg_{g}_AF"] is not None]
    flags = []
    if all(v is None for v in af.values()):
        flags.append("NO_FREQ")
    elif grpmax is None:
        flags.append("NO_GRPMAX")
    if excl and max(excl) >= MAX_AF:
        flags.append("COMMON_IN_EXCLUDED_GROUP")
    if af["gnomADg_AF"] is not None and af["gnomADg_AF"] >= MAX_AF:
        flags.append("GLOBAL_AF_COMMON")
    every_group = excl + list(grp.values())
    extra = {
        "global_rule_keeps": af["gnomADg_AF"] is None or af["gnomADg_AF"] < MAX_AF,
        "all_groups_rule_keeps": not every_group or max(every_group) < MAX_AF,
        "bub1b_1_to_5pct": "BUB1B" in symbols and grpmax is not None and MAX_AF <= grpmax < BUB1B_SENSITIVITY_AF,
        "removed_only_by_mid": pop == "MID" and grpmax >= MAX_AF and all(v < MAX_AF for g, v in grp.items() if g != "MID"),
    }
    if controls.get((cols[0], cols[1], cols[3], cols[4])):
        flags.append("POSITIVE_CONTROL")
        return "kept_positive_control", flags, grpmax, pop, extra
    if grpmax is None:
        return "kept_no_freq", flags, grpmax, pop, extra
    if grpmax >= MAX_AF:
        return "removed_common", flags, grpmax, pop, extra
    return "kept", flags, grpmax, pop, extra


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Filter 2 output VCF (filter2.vcf or filter2_no_cadd.vcf)")
    ap.add_argument("output", help="output VCF (uncompressed)")
    ap.add_argument("--controls", default=str(Path(__file__).resolve().parents[1] / "config" / "positive_controls_grouped.tsv"))
    args = ap.parse_args()

    t0 = time.time()
    controls = load_controls(Path(args.controls))
    counts = {f"{scope}_{o}": 0 for scope in ("outside", "window") for o in OUTCOMES}
    alt = {"global_rule_keeps": 0, "all_groups_rule_keeps": 0, "bub1b_1_to_5pct": 0, "removed_only_by_mid": 0}
    names = None
    n_in = 0

    with open(args.input, encoding="utf-8", errors="surrogateescape", newline="") as fin, \
         open(args.output, "w", encoding="utf-8", errors="surrogateescape", newline="") as fout:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("##INFO=<ID=CSQ,"):
                    names = line.split("Format: ", 1)[1].split('"', 1)[0].split("|")
                if line.startswith("#CHROM"):
                    if names is None:
                        sys.exit("ERROR: no ##INFO=<ID=CSQ> header")
                    fout.writelines(header_lines(Path(args.controls).name))
                fout.write(line)
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            outcome, flags, grpmax, pop, extra = classify(cols, names, controls)
            scope = "window" if "MVA_WINDOW" in cols[7] else "outside"
            counts[f"{scope}_{outcome}"] += 1
            for k, v in extra.items():
                alt[k] += v
            if outcome.startswith("removed"):
                continue
            if grpmax is not None:
                cols[7] += f";F3_GRPMAX={grpmax:g};F3_GRPMAX_POP={pop}"
            if flags:
                cols[7] += ";F3_FLAGS=" + ",".join(flags)
            fout.write("\t".join(cols) + "\n")

    kept = sum(v for k, v in counts.items() if "_kept" in k)
    summary = [f"records_in\t{n_in}", f"kept\t{kept}", f"removed\t{n_in - kept}",
               f"positive_controls_loaded\t{len(controls)}"]
    summary += [f"{k}\t{v}" for k, v in counts.items()]
    summary += [f"sensitivity_{k}\t{v}" for k, v in alt.items()]
    print("\n".join(summary), file=sys.stderr)
    print(f"elapsed_seconds\t{time.time() - t0:.0f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

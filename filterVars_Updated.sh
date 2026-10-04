set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 INPUT.vcf[.gz] OUTPUT.vcf [GENES_TSV] [WINDOW]" >&2
    exit 1
fi
INPUT=$1
OUTPUT=$2
GENES=${3:-"$(dirname "$0")/mva_genes.tsv"}
WINDOW=${4:-100000}
SECONDS=0

zcat -f "$INPUT" | awk -F'\t' -v OFS='\t' -v W="$WINDOW" '
BEGIN {
    split("splice_region_variant splice_donor_5th_base_variant splice_donor_region_variant splice_polypyrimidine_tract_variant", s, " ")
    for (k in s) SPLICE[s[k]] = 1
    nreason = split("HIGH_MODERATE SPLICE MVA_WINDOW CLINVAR_PLP CLINVAR_CONFLICTING", REASON, " ")
}

# 1st input: gene coordinate table -> windows per chromosome
FNR == NR {
    if ($0 ~ /^#/ || $0 ~ /^[ \t]*$/) next
    nw[$2]++
    lo[$2, nw[$2]] = $3 - W
    hi[$2, nw[$2]] = $4 + W
    next
}

# 2nd input: the VCF
/^#/ {
    if (index($0, "##INFO=<ID=CSQ,") == 1) {
        f = $0
        sub(/.*Format: /, "", f)
        sub(/">.*/, "", f)
        nf = split(f, names, "|")
        for (k = 1; k <= nf; k++) col[names[k]] = k
        I_IMP = col["IMPACT"]; I_CAN = col["CANONICAL"]; I_CSQ = col["Consequence"]; I_SIG = col["CLIN_SIG"]
    } else if (index($0, "#CHROM") == 1) {
        if (!I_IMP) { print "ERROR: no ##INFO=<ID=CSQ> header found" > "/dev/stderr"; exit 1 }
        print "##INFO=<ID=F1_REASON,Number=.,Type=String,Description=\"Filter 1 (annotation) rules that kept this record: HIGH_MODERATE, SPLICE, MVA_WINDOW, CLINVAR_PLP, CLINVAR_CONFLICTING\">"
        print "##Filter1=<Version=\"1.0\",Window=" W ",Spec=\"docs/07_FILTER1_AUDIT.md section E\",CSQ=\"unchanged\",QualityRarity=\"not applied\">"
    }
    print
    next
}

{
    n_in++
    pos = $2 + 0
    inwin = 0
    for (k = 1; k <= nw[$1]; k++) if (pos >= lo[$1, k] && pos <= hi[$1, k]) { inwin = 1; break }

    # Fast path: every annotation rule needs one of these words somewhere in the line.
    if (!inwin && $0 !~ /HIGH|MODERATE|splice_|pathogenic/) next

    hm = sp = plp = conf = 0
    csq = ""
    ninfo = split($8, info, ";")
    for (k = 1; k <= ninfo; k++) if (index(info[k], "CSQ=") == 1) { csq = substr(info[k], 5); break }
    if (csq != "") {
        nann = split(csq, ann, ",")
        for (j = 1; j <= nann; j++) {
            split(ann[j], a, "|")
            if (a[I_CAN] == "YES") {
                if (a[I_IMP] == "HIGH" || a[I_IMP] == "MODERATE") hm = 1
                if (!sp) {
                    nt = split(a[I_CSQ], terms, "&")
                    for (t = 1; t <= nt; t++) if (terms[t] in SPLICE) { sp = 1; break }
                }
            }
            if (a[I_SIG] != "") {
                nwd = split(a[I_SIG], words, /[&\/]/)
                for (t = 1; t <= nwd; t++) {
                    if (words[t] == "pathogenic" || words[t] == "likely_pathogenic") plp = 1
                    else if (words[t] == "conflicting_classifications_of_pathogenicity") conf = 1
                }
            }
        }
    }

    hit[1] = hm; hit[2] = sp; hit[3] = inwin; hit[4] = plp; hit[5] = conf
    reasons = ""
    for (r = 1; r <= nreason; r++) if (hit[r]) { reasons = reasons (reasons == "" ? "" : ",") REASON[r]; cnt[r]++ }
    if (reasons == "") next

    $8 = $8 ";F1_REASON=" reasons
    print
    n_out++
}

END {
    printf "records_in\t%d\nrecords_kept\t%d\n", n_in, n_out > "/dev/stderr"
    for (r = 1; r <= nreason; r++) printf "reason_%s\t%d\n", REASON[r], cnt[r] > "/dev/stderr"
}
' "$GENES" - > "$OUTPUT"

echo -e "elapsed_seconds\t$SECONDS" >&2

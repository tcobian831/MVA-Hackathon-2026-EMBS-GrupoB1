#!/bin/bash

#First filter the variants with a HIGH and MODERATE impact classification.
grep -E 'HIGH|MODERATE' vep/WGS_EX2312012_HGWCNDSX7.completo.vep.vcf > HiMo_imp_var.vep

# You should get 19353 genomic positions.
cat HiMo_imp_var.vep | wc -l

#For each position there are 1 to 600 transcripts, so we need to filter for status 'Canonical transcript'. Note that the vep file has columns separated by '\t', in the 8th column you could find general info for the genomic position organized into columns separated by a ';', at the last column you  could find the VEP anotation (CSQ=), but in there are differents transcrips separated by ','. The information of each transcript is then organizated in columns separated by '|', the column 24 corresponds to the canonical status for tha transcript. We filter only the transcripst that have the label 'YES' in that column, using awk:
awk -F'\t' '
BEGIN {
    OFS = "\t"
}
!/^#/ {

    n_info = split($8, info, ";");
    nueva_info = "";


    for (i = 1; i <= n_info; i++) {

        if (info[i] ~ /^CSQ=/) {

            sub(/^CSQ=/, "", info[i]);


            n_trans = split(info[i], trans, ",");
            nueva_csq = "CSQ=";
            first_kept = 1;


            for (j = 1; j <= n_trans; j++) {

                m = split(trans[j], campos, "|");


                if (m >= 24 && campos[24] == "YES") {
                    if (first_kept) {
                        nueva_csq = nueva_csq trans[j];
                        first_kept = 0;
                    } else {
                        nueva_csq = nueva_csq "," trans[j];
                    }
                }
            }

            info[i] = nueva_csq;
        }


        if (i == 1) {
            nueva_info = info[i];
        } else {
            nueva_info = nueva_info ";" info[i];
        }
    }

    $8 = nueva_info;
    print $0;

}' HiMo_imp_var.vep > HiMo_imp_var_Can.vep

#Then we  have to refilter only the transcripts with a HIGH and MODERATE impact classification.

awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {
    if ($8 ~ /CSQ=/) {
        n_info = split($8, info, ";");
        for (i = 1; i <= n_info; i++) {
            if (info[i] ~ /^CSQ=/) {
                sub(/^CSQ=/, "", info[i]);
                n_trans = split(info[i], trans, ",");
                nueva_csq = "CSQ=";
                first_kept = 1;
                for (j = 1; j <= n_trans; j++) {
                    m = split(trans[j], campos, "|");
                    if (m >= 3 && (campos[3] == "HIGH" || campos[3] == "MODERATE")) {
                        if (first_kept) {
                            nueva_csq = nueva_csq trans[j];
                            first_kept = 0;
                        } else {
                            nueva_csq = nueva_csq "," trans[j];
                        }
                    }
                }
                info[i] = nueva_csq;
            }
            if (i == 1) { nueva_info = info[i]; } else { nueva_info = nueva_info ";" info[i]; }
        }
        $8 = nueva_info;
    }
    print $0;
}' HiMo_imp_var_Can.vep > transCan_HiMo_imp_vars.vep


#Finally I add an extra step. I filter for low deep and genotype quality.
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ { n_fmt = split($9, fmt, ":"); n_val = split($10, val, ":"); dp = 0; gq = 0; for (i = 1; i <= n_fmt; i++) { if (fmt[i] == "DP") dp = val[i]; if (fmt[i] == "GQ") gq = val[i]; } if (dp >= 10 && gq >= 20) { print $0; } }' transCan_HiMo_imp_vars.vep > transCan_HiMo_imp_vars_qualityFilter.vep

#At the final step we should have 19011 genomic positions each one, with 2 canonical transcripts max.

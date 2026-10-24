# Additional-variant review and limits

The full VEP was reopened after ranking. The original caller VCF had 5,012,204 records; the VEP input has 4,943,948. In the earlier independent raw-VCF audit, all 2,401 records within the 15 target genes and immediate flanks were represented in the VEP comparison (no missing target-interval records). This does not claim that every original contig or every potential structural variant is represented.

## Accepted additions

| Addition | Original call | Evidence | Decision |
|---|---|---|---|
| LRP4 LR2, chr11:46878979 T>C, p.Asn1022Asp | GT 0/1; DP43/GQ99; AD23,20 | CADD19.53; groupmax0.1661%; missense, uncertain/likely benign annotation snapshot | Included by participant instruction as an exploratory partner to LR1; unknown phase. |
| TACC2 T18, chr10:122194203 CA>C | GT0/1; DP44/GQ99; AD26,18 | CADD unavailable; groupmax0.00294%; intronic deletion | Included by participant instruction as an exploratory partner to T8; no proven splice effect; unknown phase to T8. |

These are the only two variants retained as additions in the later review. They are stored separately in config/post_ranking_additions.csv and are NOT inputs to the group ranking. The earlier experimental re-ranking with them is not the reported group consensus.

## Why other additional variants were not included

The per-allele audit lists objective failures separately from manual non-inclusion. A variant passing the loose review cutoffs is not necessarily a credible second damaging allele. Non-inclusion does not establish benignity.

* BUB1B B2 and B7 are intronic/flanking, with groupmax about3.061% and CADD5.521/4.025. They were not preferred over the existing stop/missense hypothesis. No regulatory mechanism was established.
* TACC2 T2/T3 are alternate intronic repeat insertions at the same site, called1/2. One is opposite the other at that site, but their phase relative to T8 is unknown. No damaging transcript effect was established. T4 is an intronic repeat insertion with groupmax3.548%; T21 is an intronic/flanking deletion at exactly5% (eligible under <=, unlike the earlier strict< search). Neither was selected as the additional partner. T18 was the explicit chosen exploratory addition, not an experimentally validated superior allele.
* SVIL repeat insertions SV3/SV4/SV5 have frequencies around3–5% and no demonstrated damaging mechanism. SV16 is homozygous, DP7/GQ20, with AD0,4 and missing CADD; it passes the relaxed DP/GQ rule but has limited read support and no demonstrated intronic effect. SV10/SV11 are intronic/flanking with CADD5.317/6.758 and groupmax3.851%/1.004%. SV1–SV11 phase is unknown. No additional SVIL variant was approved.
* DIDO1 D2 is a homozygous intronic repeat insertion, groupmax3.798%, missing CADD. It occurs on both copies and is not a complementary second heterozygous allele to D1. No damaging mechanism was established.
* Other LZTR1 and PRMT6 additions in the earlier 88-allele search have CADD below4. CTU2 CT1 was resolved to CADD3.590 and is excluded, as are the other low-CADD CTU2 alternatives. CT4 already has GT1/1.
* TRDMT1 TR2/TR3 have CADD1.117/2.798 and fail CADD>=4. TR1 remains the only ranking variant.
* MAD2L1BP has no gene-annotated allele passing all latest criteria in the audited data. The rare intronic candidate has CADD1.017. Both historical insertion assignments were window-only, with original GT1/2 but AD7,1,0. Their local VEP target is a noncoding transcript (ENSG00000236961), not MAD2L1BP. They were removed entirely from ranking.

Read `results/expanded_review.csv` for every scanned called ALT with quality, original GT/AD/phase metadata, frozen CADD and frequency. The search covers the historical 15-gene panel, not all 28 matrix genes. Supplemental online scores are frozen to avoid database drift; missing values mean not recovered in these sources.

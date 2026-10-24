# Repository migration

Original repository: https://github.com/tcobian831/MVA-Hackathon-2026-EMBS-GrupoB1

| Original item | Disposition |
|---|---|
| filter_vars.py | 01_annotation_filter.py |
| filterVars.sh | removed: superseded initial shell experiment |
| filterVars_Updated.sh | removed: shell duplicate of filter_vars.py |
| qfilter_script.py | historical Lucia branch only; retained as documented provenance, not executed in the consolidated ranking |
| filter_population_frequency_ANY_v3.py | removed: global AF ANY alternative, not groupmax |
| final_filter/filter_final_candidates.py | 06_shortlist_filter.py |
| wmatrix.py / wmatrix_LQMV.py | replaced by matrix_engine.py plus final stage 08 |

The downloaded historical source files remain outside the deliverable repository. Drive and GitHub were not modified. This ZIP is a replacement working tree, not an additive overlay: remove superseded tracked files when adopting it.

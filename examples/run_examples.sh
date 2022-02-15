#!/bin/bash

refpanel="example_data/EUR.1KG.GRCh37.chr11.subset"
mqtl="example_data/mQTL_cg07689907"
eqtl="example_data/eQTL_FADS2"

# cg07689907 -> FADS2 -> Triglyceride example
./smrivw --mediation-cis --bfile $refpanel --gwas-summary example_data/TRI_UKBB_gwas_summary.subset.ma --beqtl-summary $mqtl --beqtl-summary $eqtl --out mqtl_eqtl_TRI_mediation --peqtl-smr 1e-06 --pmqtl-smr 1e-06 --ld-multi-snp 0.05 --thread-num 1 --diff-freq-prop 0.4 --diff-freq 0.1 --smr-wind 500 --ld-matrix --p-expo-med 0.01 --get-snp-effects

# cg07689907 -> FADS2 -> Apolipoprotein A example
./smrivw --mediation-cis --bfile $refpanel --gwas-summary example_data/APOA_UKBB_gwas_summary.subset.ma --beqtl-summary $mqtl --beqtl-summary $eqtl --out mqtl_eqtl_APOA_mediation --peqtl-smr 1e-06 --pmqtl-smr 1e-06 --ld-multi-snp 0.05 --thread-num 1 --diff-freq-prop 0.4 --diff-freq 0.1 --smr-wind 500 --ld-matrix --p-expo-med 0.01 --get-snp-effects




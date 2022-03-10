## Analysis: ##

``R CMD BATCH R_commands_proliferation_under_drug``

## Selection of analysis schemes meeting applicability criteria (residual normality and variance homogeneity): ##

``for f in `ls Applicability_*.txt`;do if test `grep -v '"' $f | awk '$2>=0.05 {print $2}' | wc -l` -eq 2;then echo $f;fi;done``

shows that only "one variance per clone" analyses meet the applicability criteria (only for assays #1 and 3). Because assay #3 appears more reliable (cells were seeded more homogeneously), Figure 3 D and E uses data from "assay 3" only (raw and processed data from assays #1 and 2 are in 'Additional_assays_not_shown.tar.bz2').

Outputs of the analysis of assay #3 are in 'Outputs.tar.bz2' (results of the analysis of variance comparing naïve and genotype-informed models are in 'Model_comparison_output_DOX_assay_3.txt' and 'Model_comparison_output_5FU_assay_3.txt' for doxorubicin and 5-fluoro-uracil respectively).

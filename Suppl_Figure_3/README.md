## Stratification of cases by tumor grade: ##

Extraction of tumor grade from file 'clinical.tsv', downloaded as in https://github.com/HKeyHKey/Mockly_et_al_2021/tree/main/Figure_1/Panel_A:

``awk -F '\t' '{print $2,$125}' clinical.tsv > Tumor_grades.dat``

(*N.B.*: some case ID's appear twice in 'clinical.tsv', but always with the same tumor grade in the two occurences)

Extraction of case ID's whose tumor grade could be assessed, and annotation of the miRNA expression data file (file 'sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv', generated in https://github.com/HKeyHKey/Mockly_et_al_2021/tree/main/Figure_1/Panel_A) with tumor grade information:

``tail -n +2 Tumor_grades.dat | awk '{print $2}' | sort | uniq -c;head -1 Tumor_grades.dat > Available_tumor_grades.dat;grep ' G[1-4]$' Tumor_grades.dat | sort | uniq >> Available_tumor_grades.dat;./Module_extracts_expression_data_with_available_tumor_grade.pl Available_tumor_grades.dat sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv > miRNA_expression_data_with_tumor_grade.tsv``

## Plot tracing: ##

``R CMD BATCH R_commands_miRNA_expression_primary_tumor_vs_NAT_by_tumor_grade``

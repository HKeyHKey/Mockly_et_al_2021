## Stratification of cases by tumor grade: ##

Extraction of tumor grade from file 'clinical.tsv', downloaded as in https://github.com/HKeyHKey/Mockly_et_al_2021/tree/main/Figure_1/Panel_A:

``awk -F '\t' '{print $2,$125}' clinical.tsv > Tumor_grades.dat``

(*N.B.*: some case ID's appear twice in 'clinical.tsv', but always with the same tumor grade in the two occurences)

Extraction of case ID's whose tumor grade could be assessed, and annotation of the miRNA expression data file (file 'sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv', generated in https://github.com/HKeyHKey/Mockly_et_al_2021/tree/main/Figure_1/Panel_A) with tumor grade information:

## Data download: ##

List of cases analyzed in Figure 1A (starting from file 'sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv', generated in https://github.com/HKeyHKey/Mockly_et_al_2021/tree/main/Figure_1/Panel_A):

``tail -n +2 sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv | awk -F '\t' '{print $3}' | sort | uniq > Cases_with_miRNA_expression_available.txt``

Selection of those cases on [the GDC portal](https://portal.gdc.cancer.gov/exploration) (using the "Upload Case Set" option), then sub-selection of cases where the mutation status of p53 is available ([Direct link to query](https://portal.gdc.cancer.gov/exploration?filters=%7B%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.available_variation_data%22%2C%22value%22%3A%5B%22ssm%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22set_id%3ADqcEp3sBQZZJsd_oiGiD%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22ENSG00000141510%22%5D%7D%7D%5D%2C%22op%22%3A%22and%22%7D&searchTableTab=cases)) (n=1,104 cases as of August 12, 2021). Case set downloaded as 'case_set_ssm__input_set__ENSG00000141510.2021-08-12.tsv' (see 'Screenshot_\*_p53_status.png' for step-by-step description of procedure).

Case stratification by p53 mutational status (using case nomenclature from file 'clinical.tsv', downloaded as in https://github.com/HKeyHKey/Mockly_et_al_2021/tree/main/Figure_1/Panel_A):

``for id in `tail -n +2 case_set_ssm__input_set__ENSG00000141510.2021-08-12.tsv`;do awk '$1=="'$id'" {print $2}' clinical.tsv;done | sort | uniq > Cases_with_p53_mutation.txt;cat Cases_with_miRNA_expression_available.txt Cases_with_p53_mutation.txt Cases_with_p53_mutation.txt | sort | uniq -c | awk '$1==1 {print $2}' > Cases_without_p53_mutation.txt``

## Extraction of miRNA expression data files for both types of cases (intact or mutant p53): ##

``./Module_extracts_expression_data_with_or_without_p53_mutation.pl Cases_with_p53_mutation.txt sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv > miRNA_expression_data_with_p53_mutation.tsv;./Module_extracts_expression_data_with_or_without_p53_mutation.pl Cases_without_p53_mutation.txt sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv > miRNA_expression_data_without_p53_mutation.tsv``

## Tracing heatmap for miRNA expression in p53-intact cases: ##

``Rscript R_commands_miRNA_expression_primary_tumor_vs_NAT_stratification_by_p53_status without``

Resulting files: heatmap in 'Heatmap_miRNA_expression_without_p53_mutation.pdf', numerical data for heatmap in 'Data_for_miRNA_expression_heatmap_without_p53_mutation.csv' and 'p-val_for_miRNA_expression_heatmap_without_p53_mutation.csv'.

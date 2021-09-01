## Data download: ##

Data downloaded from https://portal.gdc.cancer.gov/exploration on April 29, 2021, selecting "Solid tissue normal" samples in the "Exploration" tab ([Direct link to Query](https://portal.gdc.cancer.gov/exploration?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22solid%20tissue%20normal%22%5D%7D%7D%5D%7D)), then from the "Repository" tab, sub-selecting cases with "miRNA-Seq" as an experimental strategy, and "open" as an access option. Metadata related to that dataset (Biospecimen, Clinical, Sample Sheet, Metadata, and manifest) was downloaded from the GDC portal cart. Manifest for that dataset: 'gdc_manifest_20210429_152630.txt'.

Download:

``./gdc-client download -m gdc_manifest_20210429_152630.txt``

Clinical annotation un-archiving:

``tar -xzf clinical.cart.2021-04-29.tar.gz``

Data download from miRBase v.21 (which is the miRBase version used by GDC as of April 29, 2021):

``wget ftp://mirbase.org/pub/mirbase/21/high_conf_hairpin.fa.gz;gunzip high_conf_hairpin.fa.gz;mv high_conf_hairpin.fa high_conf_hairpin_miRBase21.fa``


## Data processing: ##

Extracting small RNA expression from every hairpin precursor, and annotating it with miRBase "confidence" level of the hairpin, and case and sample ID and description:

``./Module_sorted_miRNA_expression_in_primary_tumor_vs_NAT.pl gdc_sample_sheet.2021-04-29.tsv clinical.tsv``

Resulting file: 'sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv' (available here in bzip2-compressed form).

## Heatmap tracing: ##

``R CMD BATCH R_commands_miRNA_expression_primary_tumor_vs_NAT``

Resulting files: 'Heatmap_miRNA_expression.pdf' for the graphical heatmap, 'Data_for_miRNA_expression_heatmap.csv' and 'p-val_for_miRNA_expression_heatmap.csv' for its numeric data.

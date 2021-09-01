## Data download: ##

Data downloaded from https://portal.gdc.cancer.gov/exploration on February 24, 2021, selecting "maf" as data format and "open" access ([Direct link to query](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22maf%22%5D%7D%7D%5D%7D)). Manifest downloaded as 'gdc_manifest.2021-02-24.txt' and clinical data and sample sheet (downloaded from the GDC portal cart) saved as 'clinical.tsv' and 'gdc_sample_sheet.2021-02-24.tsv' respectively.

Download:

``gdc-client download -m gdc_manifest.2021-02-24.txt``

Data download from [miRBase](https://www.mirbase.org/) v.21 (which is the miRBase version used by GDC as of April 29, 2021):

``wget ftp://mirbase.org/pub/mirbase/21/high_conf_hairpin.fa.gz;wget ftp://mirbase.org/pub/mirbase/21/high_conf_mature.fa.gz;gunzip high_conf_hairpin.fa.gz high_conf_mature.fa.gz;mv high_conf_hairpin.fa hairpin_high_conf.fa;mv high_conf_mature.fa mature_high_conf.fa;wget ftp://mirbase.org/pub/mirbase/21/genomes/hsa.gff3``

## Data processing: ##

``./Module_GDC_cancer_simple_nucleotide_variation_in_miRNAs.pl clinical.tsv gdc_sample_sheet.2021-02-24.tsv 100``

Resulting file: 'Observed_cancer_simple_nucleotide_variation_in_hsa_miRNAs_cutoff_100.tsv' (available here in bzip2-compressed form).

## Heatmap tracing: ##

``Rscript R_commands_heatmap_cancer_variation_density_in_miRNAs Observed_cancer_simple_nucleotide_variation_in_hsa_miRNAs_cutoff_100.tsv``

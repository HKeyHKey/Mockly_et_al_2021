## Data download: ##

Data downloaded from https://portal.gdc.cancer.gov/exploration on February 24, 2021, selecting "maf" as data format and "open" access ([Direct link to query](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22maf%22%5D%7D%7D%5D%7D)). Manifest downloaded as 'gdc_manifest.2021-02-24.txt' and clinical data and sample sheet (downloaded from the GDC portal cart) saved as 'clinical.tsv' and 'gdc_sample_sheet.2021-02-24.tsv' respectively.

Download:

``gdc-client download -m gdc_manifest.2021-02-24.txt``

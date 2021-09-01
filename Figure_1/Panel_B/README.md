## Data download: ##

Data downloaded from https://portal.gdc.cancer.gov/repository on March 4, 2021, selecting Data Category as "copy number variation", Data Type as "Gene Level Copy Number" and Access as "open" ([Direct link to query](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22copy%20number%20variation%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Level%20Copy%20Number%22%5D%7D%7D%5D%7D)):

``mkdir Unedited;cd Unedited;gdc-client download -m ../gdc_manifest.2021-03-04.txt``

Metadata download from the GDC portal cart: clinical data and Sample sheeet saved as 'clinical.tsv' and 'gdc_sample_sheet.2021-03-04.tsv' respectively.

High-confidence miRNA hairpin data downloaded from [miRBase](https://www.mirbase.org/) v.22.1:

``wget ftp://mirbase.org/pub/mirbase/22.1/hairpin_high_conf.fa.gz;gunzip hairpin_high_conf.fa.gz``


## Data correction: ##

In the GDC copy number files, 16 miRNA hairpins are annotated at several distinct genomic loci (but with the same hairpin name); *e.g.*, MIR3670-2 appears at chr16:16306370-16306434 (the real locus for mir-3670-2) and at chr16:18488301-18488365 (whereas this is the locus for mir-3670-4). Manual correction of the GDC files for these 16 hairpins (using [miRBase](https://www.mirbase.org/) v.22 as a reference):

``cd ..;wget ftp://mirbase.org/pub/mirbase/22/genomes/hsa.gff3;for file in `ls Unedited/*/*.tsv`;do dest=`echo $file | sed 's|^Unedited/||'`;dir=`echo $dest | sed 's|/.*||'`;mkdir $dir;sed -e '/^ENSG00000258354.1\tMIR3180-1\tchr16\t14909887\t14911345\t/ d' -e '/^ENSG00000265373.1\tMIR3180-1\tchr16\t16309875\t16309968\t/ s|MIR3180-1|MIR3180-2|' -e '/^ENSG00000257366.1\tMIR3180-2\tchr16\t16308542\t16310000\t/ d' -e '/^ENSG00000266291.1\tMIR3180-2\tchr16\t18402178\t18402271\t/ s|MIR3180-2|MIR3180-3|' -e '/^ENSG00000265537.1\tMIR3180-3\tchr16\t14911220\t14911313\t/ s|MIR3180-3|MIR3180-1|' -e '/^ENSG00000257563.1\tMIR3180-3\tchr16\t18402146\t18403604\t/ d' -e '/^ENSG00000257391.1\tMIR3180-4\tchr16\t15154903\t15157020\t/ d' -e '/^ENSG00000275391.1\tMIR6770-2\tchr16\t18379351\t18379410\t/ s|MIR6770-2|MIR6770-3|' -e '/^ENSG00000277698.1\tMIR6770-1\tchr16\t16329305\t16329364\t/ s|MIR6770-1|MIR6770-2|' -e '/^ENSG00000265658.4\tMIR3690\tchrX\t1293918\t1293992\t/ s|MIR3690|MIR3690-1|' -e '/^ENSGR0000265658.4\tMIR3690\tchrY\t1293918\t1293992\t/ s|MIR3690|MIR3690-2|' -e '/^ENSG00000277120.3\tMIR6089\tchrX\t2609191\t2609254\t/ s|MIR6089|MIR6089-1|' -e '/^ENSGR0000277120.3\tMIR6089\tchrY\t2609191\t2609254\t/ s|MIR6089|MIR6089-2|' -e '/^ENSG00000274301.1\tMIR3179-3\tchr16\t14901508\t14901591\t/ s|MIR3179-3|MIR3179-1|' -e '/^ENSG00000257381.3\tMIR3179-4\tchr16\t16300159\t16300242\t/ s|MIR3179-4|MIR3179-2|' -e '/^ENSG00000274301.1\tMIR3179-3\tchr16\t14901508\t14901591\t/ d' -e '/^ENSG00000257527.1\tMIR3179-3\tchr16\t18411309\t18411851\t/ d' -e '/^ENSG00000277014.1\tMIR3179-1\tchr16\t18494493\t18494576\t/ s|MIR3179-1|MIR3179-4|' -e '/^ENSG00000257264.3\tMIR3179-1\tchr16\t14901499\t14902174\t/ d' -e '/^ENSG00000275105.1\tMIR6511B1\tchr16\t15134066\t15134150\t/ d' -e '/^ENSG00000275259.1\tMIR6511A2\tchr16\t14925937\t14926003\t/ d' -e '/^ENSG00000276311.1\tMIR6511A2\tchr16\t16368876\t16368942\t/ d' -e '/^ENSG00000272920.1\tMIR1253\tchr17\t2748078\t2748182\t/ d' -e '/^ENSG00000265496.3\tMIR1539\tchr18\t49487339\t49491878\t/ d' -e '/^ENSG00000274025.1\tMIR3670-2\tchr16\t18488301\t18488365\t/ s|MIR3670-2|MIR3670-4|' -e '/^ENSG00000221190.1\tMIR1184-1\tchrX\t155383100\t155383198\t/ s|MIR1184-1|MIR1184-2|' -e '/^ENSG00000275950.1\tMIR6724-3\tchr21\t8205315\t8205406\t/ s|MIR6724-3|MIR6724-1|' -e '/^ENSG00000275692.1\tMIR6724-3\tchr21\t8432530\t8432621\t/ s|MIR6724-3|MIR6724-4|' $file > $dest;done``

## Data processing: ##

``./Module_plots_copy_number_per_cancer_type.pl gdc_sample_sheet.2021-03-04.tsv clinical.tsv``

## Heatmap tracing: ##

``R CMD BATCH R_commands_heatmap_cancer_copy_number_in_miRNAs``

Processing table with numeric values of the heatmap:

``head -1 Data_for_copy_number_heatmap.csv > tmp_header.csv;tail -n +2 Data_for_copy_number_heatmap.csv | awk -F ',' '{print $1}' | sed -e 's|\.\.NOS|, NOS|g' -e 's|\.\.\([A-Z]\)| (\1|' -e 's|\."|)"|' -e 's|\.\.|, |g' -e 's|\.| |g' > tmp_column1.csv;tail -n +2 Data_for_copy_number_heatmap.csv | sed 's|.*",||' > tmp_other_columns.csv;mv tmp_header.csv Data_for_copy_number_heatmap.csv;paste -d ',' tmp_column1.csv tmp_other_columns.csv >> Data_for_copy_number_heatmap.csv``

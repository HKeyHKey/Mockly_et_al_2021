## Data download: ##

Small RNA-Seq described in [Isakova et al., 2020](https://www.pnas.org/doi/abs/10.1073/pnas.2002277117?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed): SRA accession numbers listed in file 'List_SRX.txt'.

Data files downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra) using fastq-dump.

## Extraction of miRNA read counts and sequencing depths: ##

``./Script_human_tissue_Small_RNA-Seq.sh``

Resulting files: stored in archive 'miRNA_counts_in_human_tissues_and_body_fluids.tar.bz2' for miRNA read counts; summarized in 'Depths.dat' for sequencing depths.

## Graph tracing: ##

``R CMD BATCH R_commands_plots_miRNA_abundance_in_tissues_bar_graphs``

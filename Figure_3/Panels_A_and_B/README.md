## Data for miRNA annotation: ##

Download of hairpin precursor sequences from miRBase v.22.1 (dated October 2018), extraction of human sequences, generation of a Bowtie2 index, and annotation of the loop and arm positions in human hairpins:

``./Script_pre-miRNA_index.sh``

Resulting files: 'hsa\_hairpinOct18.\*' (fasta and Bowtie2 index for human hairpins), 'Annotated_structures_from_hsa_hairpinOct18.dat' (annotation of the loop and arm positions).

##  miRNA quantification: ##

[Small RNA-Seq raw data](https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=695193) (in fastq format) can be downloaded from [NCBI's SRA database](https://www.ncbi.nlm.nih.gov/sra) under accession number [PRJNA695193](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA695193/) (md5sum's of the gzip-compressed files are in 'md5sum_all.txt').

Analysis of HCT-116 WT and KO clones (note that the SRA dataset also contains Small RNA-Seq data from HAP1 WT and KO clones, which are not shown in the article; they can be analyzed with the same scripts):

``for f in HCT-116KO HCT-116WT;do ./Script_miRNA_abundance_in-house_libraries.sh $f;done``

Resulting files: 'miRNA\_count\_in\_Mapping\_no\_q\_HCT-116\*.dat' files (raw number of reads for each miRNA in each genotype) and 'Depths\_no\_q.dat' (sequencing depth for each library).

## Scatter plot (Figure 3A) and barplot (Figure 3B) display: ##

``R CMD BATCH R_commands_miRNA_expression_after_mutation``

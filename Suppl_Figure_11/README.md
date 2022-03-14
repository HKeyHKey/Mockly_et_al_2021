## Data download: ##

Small RNA-Seq described in [Isakova et al., 2020](https://www.pnas.org/doi/abs/10.1073/pnas.2002277117?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed): SRA accession numbers listed in file 'List_SRX.txt'.

SRA "run" (SRR\*) accession number determination:

``./Script_download_murine_tissue_small_RNA-Seq.sh``

SRR data files downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra) using fastq-dump.

## Extraction of miRNA read counts and sequencing depths: ##

``./Script_murine_tissue_Small_RNA-Seq.sh``

Resulting files: stored in archive 'miRNA_counts_in_murine_tissues.tar.bz2' for miRNA read counts; summarized in 'Depths.dat' for sequencing depths.

## Graph tracing: ##

``sed 's|GEO accession GSM4213272 is currently private and is scheduled to be released on Dec 11, 2019.|Testis_13|' Murine_tissue_SRR_numbers.dat > Corrected_Murine_tissue_SRR_numbers.dat;for organ in `tail -n +2 Corrected_Murine_tissue_SRR_numbers.dat | awk '{print $1}' | sed 's|_[0-9]*$||' | sort | uniq`;
do Rscript R_commands_plots_miRNA_abundance_in_tissues_bar_graphs $organ `grep '^'$organ'_[0-9]* ' Corrected_Murine_tissue_SRR_numbers.dat | awk '{print $4}' | perl -pe 's/\n/ /g'`;done``

Every organ was analyzed in various biological replicates, and each biological replicate was analyzed as two technical replicates:

``for organ in `tail -n +2 Corrected_Murine_tissue_SRR_numbers.dat | awk '{print $1}' | sed 's|_[0-9]*$||' | sort | uniq`;do grep '^'$organ'_' Corrected_Murine_tissue_SRR_numbers.dat | awk '{print $1}' | sort | uniq -c;done | grep -v '^ *2 '``

Hence, because the number of technical replicates is the same for every biological replicate, averaging (without weight) every technical replicates gives the average for biological replicates.

## Localization of the loop in murine pre-miRNA hairpin structures: ##

Download of hairpin sequences, extraction of murine sequences:

``wget -O hairpinOct18.fa.gz https://www.mirbase.org/ftp/22.1/hairpin.fa.gz;gunzip hairpinOct18.fa.gz;./Fuses_lines_clean.pl hairpinOct18.fa | grep -A 1 '^>mmu\-' | grep -v '^\-\-$' > mmu_hairpinOct18.fa``

Localization of loops (relaxed parameters, trimming 15 nt on each end of the sequence before folding):

``./Module_hairpin_partitioning.pl mmu_hairpinOct18.fa 15 > mmu_problematic_hairpins.txt``

With that, 1193 hairpins were folded into unbranched hairpins and 41 were not. Folding the 41 problematic ones with a longer trimming:

``for RNA in `awk '{print $3}' mmu_problematic_hairpins.txt`;do grep -A 1 -w $RNA mmu_hairpinOct18.fa;done > missing.fa;./Module_hairpin_partitioning.pl missing_mmu.fa 25 > mmu_problematic_hairpins.txt;tail -n +2 Annotated_structures_from_missing_mmu.dat >> Annotated_structures_from_mmu_hairpinOct18.dat
``

With that, an additional 34 hairpins could be folded into unbranched hairpins. The remaining 7 had to be entered by hand, their structures are too unstable for a unique set of parameters to fold all of them properly :

``echo "mmu-mir-217 34 56 71 93" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-690 8 28 83 104" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-702 10 30 88 109" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-467f 56 76 81 101" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-1198 42 63 80 101" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-1983 12 26 100 120" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-12183 1 18 30 48" >> Annotated_structures_from_mmu_hairpinOct18.dat``


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

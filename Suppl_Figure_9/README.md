## Localization of the loop in human pre-miRNA hairpin structures: ##

Download of hairpin sequences, extraction of human sequences:

``wget -O hairpinOct18.fa.gz https://www.mirbase.org/ftp/22.1/hairpin.fa.gz;gunzip hairpinOct18.fa.gz;./Fuses_lines_clean.pl hairpinOct18.fa | grep -A 1 '^>hsa\-' | grep -v '^\-\-$' > hsa_hairpinOct18.fa``

Localization of loops (relaxed parameters, trimming 15 nt on each end of the sequence before folding):

``./Module_hairpin_partitioning.pl hsa_hairpinOct18.fa 15 > problematic_hairpins.txt``

With that, 1882 hairpins were folded into unbranched hairpins and 35 were not. Folding the 35 problematic ones with a longer trimming:

``for RNA in `awk '{print $3}' problematic_hairpins.txt`;do grep -A 1 -w $RNA hsa_hairpinOct18.fa;done > missing.fa;./Module_hairpin_partitioning.pl missing.fa 25;tail -n +2 Annotated_structures_from_missing.dat >> Annotated_structures_from_hsa_hairpinOct18.dat``

With that, an additional 29 hairpins could be folded into unbranched hairpins. The remaining 6 had to be entered by hand, their structures are too unstable for a unique set of parameters to fold all of them properly :

``echo "hsa-mir-2054 1 18 31 49" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-4536-1 1 25 64 88" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-4736 1 22 31 47" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-6134 1 57 70 109" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-6753 1 22 137 146" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-6843 1 98 102 151" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-10524 1 17 24 47" >> Annotated_structures_from_hsa_hairpinOct18.dat``

## Small RNA-Seq data download: ##

For A549 cells:

``for lib in SRR6713501 SRR5689166 DRR036695 SRR1304309;do fastq-dump $lib;done``

For Caco-2 cells:

``for lib in ERR3415707;do fastq-dump $lib;done``

For H1299 cells:

``for lib in DRR036697;do fastq-dump $lib;done``

For H2170 cells:

``for lib in SRR3341761 SRR3341762;do fastq-dump $lib;done``

For HCT-116 cells:

``for lib in SRR954987 SRR954996 SRR4235725 SRR4235726 ERR3173398 ERR3173397 ERR3173396;do fastq-dump $lib;done``

For HEK293 cells:

``for lib in SRR1240816 SRR1240817;do fastq-dump $lib;done``

For HeLa cells:

``for lib in SRR8311268 SRR8311269;do fastq-dump $lib;done``

For Hep G2 cells:

``for lib in SRR12054851;do fastq-dump $lib;done``

For IMR90 cells:

``for lib in SRR020286;do fastq-dump $lib;done``

For Kelly cells:

``for lib in SRR3533075 SRR3533074 SRR3533073;do fastq-dump $lib;done``

For RKO cells:

``for lib in ERR3415712;do fastq-dump $lib;done``

For SW480 cells:

``for lib in SRR3923807 SRR3923808;do fastq-dump $lib;done``

For U-2 OS cells:

``for lib in SRR10225092;do fastq-dump $lib;done``


## Extraction of miRNA read counts and sequencing depths: ##

``for acc in `ls *.fastq | sed 's|\.fastq$||'`;do ./Script_miRNA_abundance.sh $acc;done``

Resulting files: stored in archive 'miRNA_counts_in_human_tissues_and_body_fluids.tar.bz2' for miRNA read counts; summarized in 'Depths.dat' for sequencing depths.

## Graph tracing: ##

``for id in `ls miRNA_count_in_Mapping_SRR* | sed -e 's|miRNA_count_in_Mapping_||' -e 's|\.dat$||'`;do grep -w $id List_of_SRA_runs_for_*;done | sed -e 's|\.txt:SRR[0-9]*$||' -e 's|^List_of_SRA_runs_for_||' | sort | uniq > Human_tissues_with_Small_RNA-Seq;for tissue in `cat Human_tissues_with_Small_RNA-Seq`
do Rscript R_commands_plots_miRNA_abundance_in_tissues_bar_graphs $tissue `cat List_of_SRA_runs_for_$tissue'.txt' | perl -pe 's/\n/ /g'`
done``

## Data download: ##

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

``for acc in `ls $.fastq | sed 's|\.fastq$||'`;do ./Script_miRNA_abundance.sh $acc;done``

Resulting files: stored in archive 'miRNA_counts_in_human_tissues_and_body_fluids.tar.bz2' for miRNA read counts; summarized in 'Depths.dat' for sequencing depths.

## Graph tracing: ##

``for id in `ls miRNA_count_in_Mapping_SRR* | sed -e 's|miRNA_count_in_Mapping_||' -e 's|\.dat$||'`;do grep -w $id List_of_SRA_runs_for_*;done | sed -e 's|\.txt:SRR[0-9]*$||' -e 's|^List_of_SRA_runs_for_||' | sort | uniq > Human_tissues_with_Small_RNA-Seq;for tissue in `cat Human_tissues_with_Small_RNA-Seq`
do Rscript R_commands_plots_miRNA_abundance_in_tissues_bar_graphs $tissue `cat List_of_SRA_runs_for_$tissue'.txt' | perl -pe 's/\n/ /g'`
done``

## Data download: ##

``./Script_download_human_tissue_small_RNA-Seq.sh``

Resulting file: 'Small_RNA-Seq_in_healthy_human_tissues.tsv'.

## Data selection (because 'Small_RNA-Seq_in_healthy_human_tissues.tsv' still contains weird things: HeLa cells, "not collected" samples, ...): ##

``for tissue in 'Peripheral blood' Bone BREAST Colon LUNG 'not collected' 'Peripheral blood mononuclear cell' PROSTATE saliva Semen Thyroid 'total WBC';do match=`echo $tissue | sed 's| |_|g'`;   display=$match;case "$tissue" in "BREAST") display='Breast';;"LUNG") display='Lung';;"PROSTATE") display='Prostate';;"total WBC") display='Total_white_blood_cells';;"not collected") display='Astrocytes';;esac
echo "display="$display;sed 's| |_|g' Small_RNA-Seq_in_healthy_human_tissues.tsv | awk -F '\t' '$3=="'$match'" {print $2}' > List_of_SRA_runs_for_$display'.txt';done``

## Extraction of miRNA read counts and sequencing depths: ##

After download of the 'SRR\*.fastq' files from [SRA](https://www.ncbi.nlm.nih.gov/sra)

Resulting files: stored in archive 'miRNA_counts_in_human_tissues_and_body_fluids.tar.bz2' for miRNA read counts; summarized in 'Depths.dat' for sequencing depths.

## Graph tracing: ##

``R CMD BATCH R_commands_plots_miRNA_abundance_in_tissues_bar_graphs``

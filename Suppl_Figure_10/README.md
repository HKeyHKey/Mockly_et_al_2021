## Data download: ##

``./Script_download_human_tissue_small_RNA-Seq.sh``

Resulting file: 'Small_RNA-Seq_in_healthy_human_tissues.tsv'.

## Data selection (because 'Small_RNA-Seq_in_healthy_human_tissues.tsv' still contains weird things: HeLa cells, "not collected" samples, ...): ##

``for tissue in 'Peripheral blood' Bone BREAST Colon LUNG 'not collected' 'Peripheral blood mononuclear cell' PROSTATE saliva Semen Thyroid 'total WBC';do match=`echo $tissue | sed 's| |_|g'`;   display=$match;case "$tissue" in "BREAST") display='Breast';;"LUNG") display='Lung';;"PROSTATE") display='Prostate';;"total WBC") display='Total_white_blood_cells';;"not collected") display='Astrocytes';;esac;echo "display="$display;sed 's| |_|g' Small_RNA-Seq_in_healthy_human_tissues.tsv | awk -F '\t' '$3=="'$match'" {print $2}' > List_of_SRA_runs_for_$display'.txt';done``

## Extraction of miRNA read counts and sequencing depths: ##

After download of the 'SRR\*.fastq' files from [SRA](https://www.ncbi.nlm.nih.gov/sra) using fastq-dump:

``./Script_human_tissue_Small_RNA-Seq.sh``

(N.B.: that script uses file 'Annotated_structures_from_hsa_hairpinOct18.dat', generated as described in https://github.com/HKeyHKey/Mockly_et_al_2022/tree/main/Suppl_Figure_9)

Resulting files: stored in archive 'miRNA_counts_in_human_tissues_and_body_fluids.tar.bz2' for miRNA read counts; summarized in 'Depths.dat' for sequencing depths.

## Graph tracing: ##

``for id in `ls miRNA_count_in_Mapping_SRR* | sed -e 's|miRNA_count_in_Mapping_||' -e 's|\.dat$||'`;do grep -w $id List_of_SRA_runs_for_*;done | sed -e 's|\.txt:SRR[0-9]*$||' -e 's|^List_of_SRA_runs_for_||' | sort | uniq > Human_tissues_with_Small_RNA-Seq;for tissue in `cat Human_tissues_with_Small_RNA-Seq`
do Rscript R_commands_plots_miRNA_abundance_in_tissues_bar_graphs $tissue `cat List_of_SRA_runs_for_$tissue'.txt' | perl -pe 's/\n/ /g'`
done``

#!/bin/sh

species='mmu'
directory='/poolzfs/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/' # to be adapted with the location of your Bowtie2 index files
MIN_LENGTH=18
MAX_LENGTH=30
for acc in `cat ID_list_2`
do adapter=TGGAATTCTCGGGTGCCAAGG
   cutadapt -a $adapter SRA_download/$acc'.fastq' -o trimmed_$acc'.fastq' --discard-untrimmed -m $MIN_LENGTH -M $MAX_LENGTH
   bowtie2 --no-unal -x $species'_hairpinOct18' -U trimmed_$acc'.fastq' -S Mapping_$acc'.sam' > mapping_$acc'.log'
   bowtie2 --no-unal -x $directory'genome' -U trimmed_$acc'.fastq' -S Genomic_mapping_$acc'.sam' > genomic_mapping_$acc'.log'
   echo $acc `sed '1,/^@PG\t/ d' Genomic_mapping_$acc'.sam' | wc -l` >> Depths.dat
   ./Module_count_miRNA_reads.pl Mapping_$acc'.sam' Annotated_structures_from_$species'_hairpinOct18.dat'
done

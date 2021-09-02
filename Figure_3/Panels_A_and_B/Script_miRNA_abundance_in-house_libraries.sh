#!/bin/bash

lib=$1
adapter=GTTCAGAGTTCTACAGTCCGACGATCNNNN

directory='/poolzfs/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/' # Enter here the path to your Bowtie2 index for the hg38 genome
MIN_LENGTH=18
MAX_LENGTH=30

file=`ls $lib'_R2'[Rr]'everse_S'[1-4]'_L001_R2_001.fastq.gz' | sed 's|\.gz$||'`
gunzip $file'.gz'
#cutadapt -q 20 -g $adapter --discard-untrimmed -m $MIN_LENGTH -M $MAX_LENGTH $file | cutadapt -u -4 - > trimmed_$lib'.fastq' # if you want to screen reads on their quality (not done in the analyses shown in the article: the gain was negligible)
cutadapt -g $adapter --discard-untrimmed -m $MIN_LENGTH -M $MAX_LENGTH $file | cutadapt -u -4 - > trimmed_no_q_$lib'.fastq'


#bowtie2 --no-unal -x hsa_hairpinOct18 -U trimmed_$lib'.fastq' -S Mapping_$lib'.sam' > mapping_$lib'.log'
#bowtie2 --no-unal -x $directory'genome' -U trimmed_$lib'.fastq' -S Genomic_mapping_$lib'.sam' > genomic_mapping_$lib'.log'
bowtie2 --no-unal -x hsa_hairpinOct18 -U trimmed_no_q_$lib'.fastq' -S Mapping_no_q_$lib'.sam' > mapping_no_q_$lib'.log'
bowtie2 --no-unal -x $directory'genome' -U trimmed_no_q_$lib'.fastq' -S Genomic_mapping_no_q_$lib'.sam' > genomic_mapping_no_q_$lib'.log'

#echo $lib `sed '1,/^@PG\t/ d' Genomic_mapping_$lib'.sam' | wc -l` >> Depths.dat
echo $lib `sed '1,/^@PG\t/ d' Genomic_mapping_no_q_$lib'.sam' | wc -l` >> Depths_no_q.dat

#./Module_count_miRNA_reads.pl Mapping_$lib'.sam' Annotated_structures_from_hsa_hairpinOct18.dat
./Module_count_miRNA_reads.pl Mapping_no_q_$lib'.sam' Annotated_structures_from_hsa_hairpinOct18.dat


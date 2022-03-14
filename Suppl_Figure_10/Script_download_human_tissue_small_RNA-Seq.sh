#!/bin/sh

counter=0
retstart=0
n=100000
while test $n -eq 100000
do wget -O esearch_out_SRA_$counter https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term=Homo%20sapiens%5borgn%5d\(\"Small%20RNA-Seq\"%5ball%20fields%5b%20OR%20\"miRNA-Seq\"%5ball%20fields%5d%20OR%20mirna_seq%5bstrategy%5d%20OR%20ncrna_seq%5bstrategy%5d\)%20public%5baccess%5d%20single%5blayout%5d%20\"filetype%20fastq\"%5bfilter%5d%20transcriptomic%5bsource%5d\&retmax=100000\&retstart=$retstart
   retstart=`echo $retstart+100000 | bc`
   n=`grep -Pc '^\t*<Id>' esearch_out_SRA_$counter`
   counter=`echo $counter"+1" | bc`
   sleep 1
done

counter=0
retstart=0
n=100000
while test $n -eq 100000
do wget -O esearch_out_BioSample_$counter https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample\&term=Homo%20sapiens%5borgn%5d\"biosample%20sra\"%5bfilter%5d\&retmax=100000\&retstart=$retstart
   retstart=`echo $retstart+100000 | bc`
   n=`grep -Pc '^\t*<Id>' esearch_out_BioSample_$counter`
   counter=`echo $counter"+1" | bc`
   sleep 1
done

cat esearch_out_SRA_* | grep -P '^\t*<Id>' | sed -e 's|^\t*<Id>||' -e 's|</Id>||' > ID_list_SRA
i=1
n=200
while test $n -eq 200
do j=`echo "("$i"-1)*200+1" | bc`
   tail -n +$j ID_list_SRA | head -200 > batch_$i
   list=`perl -pe 's/\n/,/g' batch_$i | sed 's|,$||'`
   wget -O efetch_out_SRA_$i https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra\&id=$list
   n=`cat batch_$i | wc -l`
   i=`echo $i"+1" | bc`
   sleep 1
done

cat esearch_out_BioSample_* | grep -P '^\t*<Id>' | sed -e 's|^\t*<Id>||' -e 's|</Id>||' > ID_list_BioSample
i=1
n=200
while test $n -eq 200
do j=`echo "("$i"-1)*200+1" | bc`
   tail -n +$j ID_list_BioSample | head -200 > batch_$i
   list=`perl -pe 's/\n/,/g' batch_$i | sed 's|,$||'`
   wget -O efetch_out_BioSample_$i https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample\&id=$list
   n=`cat batch_$i | wc -l`
   i=`echo $i"+1" | bc`
   sleep 1
done

cat efetch_out_BioSample_* | sed 's|</BioSample>|\
|g' | grep '<Attribute attribute_name="tissue" harmonized_name="tissue" display_name="tissue">' | grep '<Id db="SRA">' | grep '<Attribute attribute_name="disease" harmonized_name="disease" display_name="disease">' | sed -e 's|.*<Id db="BioSample" is_primary="1">\([A-Z0-9]*\)</Id>.*<Id db="SRA">\([A-Z0-9]*\)</Id>.*<Attribute attribute_name="tissue" harmonized_name="tissue" display_name="tissue">\([^<]*\)</Attribute>.*<Attribute attribute_name="disease" harmonized_name="disease" display_name="disease">\([^<]*\)</Attribute>.*|\1\t\2\t\3\t\4|' -e 's|.*<Id db="BioSample" is_primary="1">\([A-Z0-9]*\)</Id>.*<Id db="SRA">\([A-Z0-9]*\)</Id>.*<Attribute attribute_name="disease" harmonized_name="disease" display_name="disease">\([^<]*\)</Attribute>.*<Attribute attribute_name="tissue" harmonized_name="tissue" display_name="tissue">\([^<]*\)</Attribute>.*|\1\t\2\t\4\t\3|' > Human_tissue_samples_in_BioSample
awk -F '\t' '$4=="Control" || $4=="normal" || $4=="healthy" || $4=="Healthy" || $4=="NA" || $4=="Normal" || $4=="unaffected" || $4=="control" || $4=="None" || $4=="healthy control" || $4=="WT" || $4=="na" || $4=="normal blood" || $4=="--" || $4=="none" || $4=="Not applicable" || $4=="n/a" || $4=="NONE" || $4=="Noncancerous tissue" || $4=="Controls" || $4=="Healthy donor" || $4=="healthy donor" || $4=="Healthy non-IBD" || $4=="N/A" || $4=="Healthy Donor" || $4=="Phenotype Control" || $4=="Unaffected" || $4=="Healthy Control" || $4=="Non-demented control" || $4=="no known cardiac diseases" || $4=="CTRL" || $4=="non-diseased" || $4=="Normal non smokers" || $4=="No" || $4=="negative" || $4=="normal controls (NC)" || $4=="NORMAL" || $4=="CONTROL" || $4=="no" || $4=="no pathologic disease" || $4=="Non-malignant tissue" || $4=="control esophagus, no disease" || $4=="Normal Tissue" || $4=="None, Control" || $4=="NO" || $4=="Healthy control" || $4=="aged normal" || $4=="young normal" || $4=="UNAFFECTED" || $4=="Normal thyroid" || $4=="Normal retina" || $4=="Normal control" || $4=="none - WT" || $4=="matched_normal" || $4=="control (normal spermatogenesis)" || $4=="Normal Control" || $4=="none known" || $4=="none described" || $4=="normal melanocyte" || $4=="normal control" || $4=="none reported" || $4=="none, normal, control" || $4=="non" || $4=="No disease"  || $4=="not applicable" || $3=="liver adjacent to HCC" {print $1"\t"$2"\t"$3}' Human_tissue_samples_in_BioSample > Healthy_human_tissue_samples_in_BioSample

cat efetch_out_SRA_* > All_SRA 
awk -F '\t' '{print $2}' Healthy_human_tissue_samples_in_BioSample | sort | uniq > id_in_BioSample
sed 's|</PRIMARY_ID>|\
|g' All_SRA | sed 's|.*<PRIMARY_ID>||' | grep '^SR' | sort | uniq > id_in_SRA 
cat id_in_SRA id_in_BioSample id_in_BioSample | sort | uniq -c > tmp_unified

awk '$1==3 {print $2}' tmp_unified > Common_ids_step1
awk '$1==1 {print $2}' tmp_unified > unique_to_SRA
awk '$1==2 {print $2}' tmp_unified > unique_to_BioSample

for f in `ls efetch_out_SRA_*`;do sed 's|</EXPERIMENT_PACKAGE>|\
|g' $f | sed -e 's|.*<PRIMARY_ID>\([DES]RR\)|\1|' -e 's|</PRIMARY_ID>.*||' | grep '^[DES]RR';done > Runs_in_SRA

for id in `cat unique_to_BioSample`
do wget -O esearch_out_run_$id https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term=$id
   ID=`grep '^<Id>' esearch_out_run_$id | sed -e 's|^<Id>||' -e 's|</Id>||'`
   sleep 0.5
   wget -O efetch_out_run_$id'_'$ID https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra\&id=$ID
   sed 's|</EXPERIMENT_PACKAGE>|\
|g' efetch_out_run_$id'_'$ID | sed -e 's|.*<PRIMARY_ID>\([DES]RR\)|\1|' -e 's|</PRIMARY_ID>.*||' | grep '^[DES]RR'
   sleep 0.5
done > Runs_in_BioSample
cat Runs_in_BioSample Runs_in_BioSample Runs_in_SRA | sort | uniq -c | awk '{print $1}' | sort | uniq -c # Result: no overlap at all! I'll have to go ahead with just the 'Common_ids_step1' Small RNA-Seq runs

for id in `cat Common_ids_step1`
do wget -O esearch_out_run_$id https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term=$id
   ID_list=`grep '^<Id>' esearch_out_run_$id | sed -e 's|^<Id>||' -e 's|</Id>||'`
   sleep 0.5
   for ID in `echo $ID_list`
   do wget -O efetch_out_run_$id'_'$ID https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra\&id=$ID
      sed 's|</EXPERIMENT_PACKAGE>|\
|g' efetch_out_run_$id'_'$ID | sed -e 's|.*<PRIMARY_ID>\([DES]RR\)|\1|' -e 's|</PRIMARY_ID>.*||' | grep '^[DES]RR'
      sleep 0.5
   done
done > Runs_in_common_ids_step1
echo "BioSample SRA_run Organ" | sed -e 's| |\t|g' -e 's|_| |g' > Small_RNA-Seq_in_healthy_human_tissues.tsv
for run in `cat Runs_in_common_ids_step1`
do study=`grep -w $run efetch_out_run_* | sed -e 's|_[0-9]*:.*||' -e 's|efetch_out_run_||' | uniq`
   awk -F '\t' '$2=="'$study'" {print $1"\t\t"$3}' Healthy_human_tissue_samples_in_BioSample | sed 's|\t\t|\t'$run'\t|'
done | grep -vw SRR4228255 | grep -vw SRR4228251 | grep -vw SRR4228303 | grep -vw SRR6081962 | grep -vw SRR6081938 | grep -vw SRR6081926 | grep -vw SRR6081924 | grep -vw SRR6081925 | grep -vw SRR4228303 >> Small_RNA-Seq_in_healthy_human_tissues.tsv # Because SRR4228255, SRR4228251, SRR4228303, SRR6081962, SRR6081938, SRR4228303 are in fact RNA-Seq (not Small RNA-Seq) datasets; SRR6081926, SRR6081924, SRR6081925 are in fact from cultured cell lines

for tissue in breast lung colon spleen bone_marrow brain cervix kidney skeletal_muscle pancreas ovary stomach # a list of organs that I don't have in the initial list: let's see if I can get some with less stringent selection criteria on disease stage
do sed 's| |_|g' Human_tissue_samples_in_BioSample | awk 'BEGIN{IGNORECASE=1} $3=="'$tissue'" {print}'
done | grep -iv leukemia | grep -iv neuroblastoma | grep -iv myeloma | grep -iv carcinoma | grep -iv cancer | grep -iv tumor | awk -F '\t' '{print $2}' > Additional_to_be_checked_for_health_status
for id in `cat Additional_to_be_checked_for_health_status`
do wget -O esearch_out_run_$id https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term=$id
   ID_list=`grep '^<Id>' esearch_out_run_$id | sed -e 's|^<Id>||' -e 's|</Id>||'`
   sleep 0.5
   for ID in `echo $ID_list`
   do wget -O efetch_out_run_$id'_'$ID https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra\&id=$ID
      sed 's|</EXPERIMENT_PACKAGE>|\
|g' efetch_out_run_$id'_'$ID | sed -e 's|.*<PRIMARY_ID>\([DES]RR\)|\1|' -e 's|</PRIMARY_ID>.*||' | grep '^[DES]RR'
      sleep 0.5
   done
done > Runs_in_BioSample_step2
cat Runs_in_BioSample_step2 Runs_in_BioSample_step2 Runs_in_SRA | sort | uniq -c | awk '$1==3 {print $2}' > Runs_in_common_ids_step2
### Then, manual control (on the SRA webserver: https://www.ncbi.nlm.nih.gov/sra/) of each of these 56 NGS runs, to verify that they indeed come from healthy human tissues, and were analyzed by Small RNA-Seq. Results in the manually-edited file 'Runs_in_common_ids_step2.tsv'. Only the last 4 samples are correct ... and they were already in 'Small_RNA-Seq_in_healthy_human_tissues.tsv'

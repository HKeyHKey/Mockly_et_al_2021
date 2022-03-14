#!/bin/sh

echo "Sample SRX_accession_number SRX_UID SRR_accession_number" > Murine_tissue_SRR_numbers.dat
for acc in `cat List_SRX.txt`
do wget -O esearch_out_$acc https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term=$acc
   sleep 1
   ID=`grep '^<Id>' esearch_out_$acc | sed -e 's|^<Id>||' -e 's|</Id>$||'`
   wget -O efetch_out_$acc https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra\&id=$ID
   sleep 1
   sed 's|sample_title="|\
sample_title="|g' efetch_out_$acc | grep '^sample_title="' | sed -e 's|^sample_title="||' -e 's|".*||' | uniq > tmp_desc_$acc
   sed 's|<PRIMARY_ID>SRR|\
SRR|g' efetch_out_$acc | grep '^SRR' | sed 's|</PRIMARY_ID>.*||' > tmp_list_$acc
   sample_desc=`cat tmp_desc_$acc`
for SRR in `cat tmp_list_$acc`
   do echo $sample_desc $acc $ID $SRR >> Murine_tissue_SRR_numbers.dat
   done
done


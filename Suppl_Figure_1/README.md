## Extraction of miR-34 family member sequences in human and in mouse: ##

``wget ftp://mirbase.org/pub/mirbase/22.1/mature.fa.gz;gunzip mature.fa.gz;mv mature.fa matureOct18.fa;for species in hsa mmu;do egrep -A 1 "^>"$species"-miR-34[abc]-5p|^>"$species"-miR-449[abc] |^>"$species"-miR-449[abc]-5p" matureOct18.fa | grep -v '\-\-' > $species'_family_members.fa';done``

For some of them, the miRBase-recorded isoform is not the most abundant one (as judged by the Small RNA-Seq profile shown on miRBase), so their sequences had to be corrected by hand: that was the case for hsa-miR-34b-5p (edited into AGGCAGUGUCAUUAGCUGAUUG) and hsa-miR-449c-5p (edited into AGGCAGUGUAUUGCUAGCGGCUGU).

## Sequence alignment: ##

``clustalw hsa_family_members.fa;clustalw mmu_family_members.fa;cat hsa_family_members.aln;tac mmu_family_members.aln``

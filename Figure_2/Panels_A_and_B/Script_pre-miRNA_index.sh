#!/bins/h

# Sequence download:
if ! [ -f hairpinOct18.fa ] # Only do this if file 'hairpinOct18.fa' does not already exist
then wget ftp://mirbase.org/pub/mirbase/22.1/hairpin.fa.gz
     gunzip hairpin.fa.gz
     mv hairpin.fa hairpinOct18.fa
fi

# Bowtie2 index building:
for species in hsa # mmu
do if test `ls $species'_hairpinOct18.1.bt2' $species'_hairpinOct18.2.bt2' $species'_hairpinOct18.3.bt2' $species'_hairpinOct18.4.bt2' $species'_hairpinOct18.rev.1.bt2' $species'_hairpinOct18.rev.2.bt2' | wc -l` -ne 6 # Only do this if Bowtie2 index files do not all already exist
   then ./Fuses_lines_clean.pl hairpinOct18.fa | grep -A 1 '^>'$species'\-' | grep -v '\-\-' | sed '/^>/ !s|U|T|g' > $species'_hairpinOct18.fa'
        bowtie2-build -f $species'_hairpinOct18.fa' $species'_hairpinOct18'
   fi
done

# Locating the loop in each hairpin:
for species in hsa # mmu
do if ! [ -f Annotated_structures_from_$species'_hairpinOct18.dat' ] # Only do this if file 'Annotated_structures_from_$species"_hairpinOct18.dat"' does not already exist
   then ./Module_hairpin_partitioning.pl $species'_hairpinOct18.fa' 15 > problematic_hairpins_$species'.txt'
        for RNA in `awk '{print $3}' problematic_hairpins_$species'.txt'`
        do grep -A 1 -w $RNA $species'_hairpinOct18.fa'
        done > missing_$species'.fa'
        ./Module_hairpin_partitioning.pl missing_$species'.fa' 25
        tail -n +2 Annotated_structures_from_missing_$species'.dat' >> Annotated_structures_from_$species'_hairpinOct18.dat'
	case "$species" in "hsa") echo "hsa-mir-2054 1 18 31 49" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-4536-1 1 25 64 88" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-4736 1 22 31 47" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-6134 1 57 70 109" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-6753 1 22 137 146" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-6843 1 98 102 151" >> Annotated_structures_from_hsa_hairpinOct18.dat;echo "hsa-mir-10524 1 17 24 47" >> Annotated_structures_from_hsa_hairpinOct18.dat;; # Six human pre-miRNA hairpins could not be folded properly with these automatic scripts: they had to be entered manually.
		           "mmu") echo "mmu-mir-217 34 56 71 93" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-690 8 28 83 104" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-702 10 30 88 109" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-467f 56 76 81 101" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-1198 42 63 80 101" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-1983 12 26 100 120" >> Annotated_structures_from_mmu_hairpinOct18.dat;echo "mmu-mir-12183 1 18 30 48" >> Annotated_structures_from_mmu_hairpinOct18.dat;; # Seven murine pre-miRNA hairpins could not be folded properly with these automatic scripts: they had to be entered manually.
        esac
   fi
done

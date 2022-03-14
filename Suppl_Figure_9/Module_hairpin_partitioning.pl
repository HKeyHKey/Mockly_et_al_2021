#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script inputs (fasta file with hairpin sequences, then flank length to be trimmed; e.g., ./Module_hairpin_partitioning.pl hsa_hairpinOct18.fa 15).\n";
}
else
{
    $radical=$ARGV[0];
    $FLANK_TRIM=$ARGV[1];
    $radical=~s/\.fa$//;
    open(OUT,">Annotated_structures_from_$radical".".dat");
    print OUT "pre-miRNA_hairpin Start_of_5p_arm End_of_5p_arm Start_of_3p_arm End_of_3p_arm\n";
    open(IN,$ARGV[0]);
    while(<IN>)
    {
        chomp;
        if (/^>/)
        {
            $name=$_;
            $name=~s/^> *//;
            $name=~s/ .*//;
        }
        else
        {
            $seq=substr $_,$FLANK_TRIM,length($_)-2*$FLANK_TRIM;
            $struct=`echo $seq | RNAfold | tail -1;rm -f rna.ps`;
            chomp $struct;
            $struct=~s/ .*//;
            $total=length($_);
            $struct=~s/\).*//;
            $arm5_loop=length($struct);
            $struct=~s/\.*$//;
            $arm5=length($struct);
            $end5=$arm5+$FLANK_TRIM;
            $start3=$arm5_loop+1+$FLANK_TRIM;
            if (($end5 < 0.3*$total) || ($start3 > 0.7*$total))
            {
                print "Problem for $name\n";
            }
            else
            {
                print OUT "$name 1 $end5 $start3 $total\n";
            }
        }
    }
    close(IN);

}

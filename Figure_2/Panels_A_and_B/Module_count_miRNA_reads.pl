#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script inputs (SAM file [mapping reads on hairpins], then text file describing hairpin secondary structures; e.g., ./Module_count_miRNA_reads.pl Mapping_SRR954987.sam Annotated_structures_from_hsa_hairpinOct18.dat).\n";
}
else
{
    open(STRUCT,$ARGV[1]);
    while(<STRUCT>)
    {
	chomp;
	if ($_ ne 'pre-miRNA_hairpin Start_of_5p_arm End_of_5p_arm Start_of_3p_arm End_of_3p_arm')
	{
		@array=split(' ',$_);
		$end_of_5p{$array[0]}=$array[2];
		$start_of_3p{$array[0]}=$array[3];
	}
    }
    close(STRUCT);

    open(IN,$ARGV[0]);
    while(<IN>)
    {
        chomp;
	if (/^\@PG\t/)
	{
		$read=1;
	}
	else
	{
		if ($read)
		{
			if (/\tMD:Z:\d+\t/)
			{
				@array=split(' ',$_);
				if ($array[1]==0)
				{
					$hairpin=$array[2];
					$position=$array[3];
					$sequence=$array[9];
					if ($position>$end_of_5p{$hairpin})
					{
						$tag='-3p';
					}
					else
					{
						if ($position+length($sequence)-1<$start_of_3p{$hairpin})
						{
							$tag='-5p';
						}
						else
						{
							$tag='-loop';
						}
					}
					++$count{$hairpin.$tag};
				}
			}
		}
	}
    }
    close(IN);
    $radical=$ARGV[0];
    $radical=~s/\.sam$//;
    open(OUT,">miRNA_count_in_$radical".".dat");
    print OUT "miRNA Number_of_reads\n";
    foreach $mature (keys %count)
    {
	print OUT "$mature $count{$mature}\n";
    }
    close(OUT);
}

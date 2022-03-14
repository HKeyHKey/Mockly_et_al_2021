#!/usr/bin/perl

if ($ARGV[0] eq '')
{
    print "Please enter input file names (e.g., ./Module_trim_5p_nucleotides_from_fastq.pl trimmed_SRR8311268.fastq).\n";
}
else
{
    $n=3;
    open(FASTQ,$ARGV[0]);
    while(<FASTQ>)
    {
	++$n;
	chomp;
	if ($n == 4)
	{
	    $name=$_;
	    $n=0;
	}
        if ($n == 1)
        {
            $l1=$_;
            $l1=~s/^....(.*)/\1/; # removes the first 4 nt
        }
        if ($n == 3)
        {
            $l3=$_;
            $l3=~s/^....(.*)/\1/; # removes the first 4 nt
            print "$name\n$l1\n+\n$l3\n";
        }
    }
    close(FASTQ);
}

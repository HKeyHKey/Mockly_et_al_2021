#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script inputs (text file with list of cases with the desired p53 mutation status, and text file with miRNA expression data; e.g., ./Module_extracts_expression_data_with_or_without_p53_mutation.pl Cases_with_p53_mutation.txt sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv).\n";
    die;
}

open(LIST,$ARGV[0]);
while(<LIST>)
{
    chomp;
    push(@select,$_);
}
close(LIST);

open(DATA,$ARGV[1]);
while(<DATA>)
{
    chomp;
    if ($_ =~ /^miRNA\tConfidence\tCase\tSample ID\tSample description\tCase description\tmiRNA rpm$/)
    {
	print "$_\n";
    }
    else
    {
	@array=split('\t',$_);
	if (grep(/^$array[2]$/,@select))
	{
	    print "$_\n";
	}
    }
}

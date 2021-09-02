#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script inputs (text file with tumor grade data, and text file with miRNA expression data; e.g., ./Module_extracts_expression_data_with_available_tumor_grade.pl Available_tumor_grades.dat sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv).\n";
    die;
}

open(LIST,$ARGV[0]);
while(<LIST>)
{
    chomp;
    if ($_ ne 'case_submitter_id tumor_grade')
    {
	@array=split(' ',$_);
	$grade{$array[0]}=$array[1];
    }
}
close(LIST);

open(DATA,$ARGV[1]);
while(<DATA>)
{
    chomp;
    if ($_ =~ /^miRNA\tConfidence\tCase\tSample ID\tSample description\tCase description\tmiRNA rpm$/)
    {
	print "$_\tTumor grade\n";
    }
    else
    {
	@array=split('\t',$_);
	if (grep(/^$array[2]$/,keys %grade))
	{
	    print "$_\t$grade{$array[2]}\n";
	}
    }
}

#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script inputs (e.g., ./Module_sorted_miRNA_expression_in_primary_tumor_vs_NAT.pl gdc_sample_sheet.2021-04-29.tsv clinical.tsv).\n";
    die;
}

open(SHEET,$ARGV[0]);
while(<SHEET>)
{
    chomp;
    if ($_ !~ '/^File ID\tFile Name\tData Category/')
    {
	@array=split('\t',$_);
	$file_path=$array[0].'/'.$array[1];
	@array2=split(', ',$array[6]);
	@array3=split(', ',$array[5]);
	@array4=split(', ',$array[7]);

	for ($n=0;$n<push(@array2);++$n)
	{
	    $sample=$array2[$n];
	    $case{$sample}=$array3[$n];
	    $sample_description{$sample}=$array4[$n];
	    push(@{$samples{$file_path}},$sample);
	}
    }
}
close(SHEET);

open(CLINICAL,$ARGV[1]);
while(<CLINICAL>)
{
    chomp;
    if ($_ !~ /^case_id\tcase_submitter_id\tproject_id\tage_at_index\t/)
    {
	@array=split('\t',$_);
	### The following command:
	# for case in `tail -n +2 clinical.tsv | awk -F '\t' '{print $2,$108,$120}' | sort | uniq | awk '{print $1}' | sort | uniq -c | grep -v '^ *1 ' | awk '{print $2}'`;do tail -n +2 clinical.tsv | awk -F '\t' '{print $2,$108,$120}' | grep -w $case;echo "";done
	### shows that every case in 'clinical.tsv' (as of Mar. 10, 2021) has a single (cancer type, cancer localization)
       	
	if (($cancer_type{$array[1]} eq '') || ($array[107] ne 'Not Reported') || ($array[107] ne 'Tumor, NOS'))
	{
	    $cancer_type{$array[1]}=$array[107];
	    $cancer_localization{$array[1]}=$array[119];
	    push(@{$cases_for_cancer_type{$array[107].' ('.$array[119].')'}},$array[1]);
	}
    }
}
close(CLINICAL);

@types=keys %cases_for_cancer_type;

for $cancer_type (@types)
{
    %seen=();
    my @unique = do { my %seen; grep { !$seen{$_}++ } @{$cases_for_cancer_type{$cancer_type}}};
    $number_of_cases{$cancer_type}=push(@unique);
}


$species='hsa';
open(CONF_HAIRPINS,'high_conf_hairpin_miRBase21.fa');
while(<CONF_HAIRPINS>)
{
    chomp;
    if (/^>$species-/)
    {
	s/^> *//;
	s/$species-//;
	s/ .*//;
	push(@high_conf_hairpins,$_);
    }
}
close(CONF_HAIRPINS);

open(MATURE,'mature_miRBase21.fa');
while(<MATURE>)
{
    chomp;
    if (/^>$species-/)
    {
	s/^> *//;
	@array=split(' ',$_);
	$human_readable_name{$array[1]}=$array[0];
    }
}
close(MATURE);


@list_files=`ls */*.isoforms.quantification.txt`; # Use "isoform" quantification rather than "miRNAs", in order to have a separate read count for miR-5p and miR-3p (in the '*.mirnaseq.mirnas.quantification.txt' files, read counts are aggregated per hairpin)

foreach $f (@list_files)
{
    chomp $f;
    push(@files,$f);
    open(DATA,$f);
#    $type_localization=$cancer_type{$case{$f}}." (".$cancer_localization{$case{$f}}.")";
    while(<DATA>)
    {
	chomp;
	if ($_ !~ /^miRNA_ID\tisoform_coords\tread_count\treads_per_million_miRNA_mapped\tcross-mapped\tmiRNA_region$/)
	{
	    @array=split('\t',$_);
	    if ($array[5]=~/^mature,/) # selects reads mapping on an arm of the hairpin
	    {
		$mimat=$array[5];
		$mimat=~s/.*,//;
		$rpm{$f}{$human_readable_name{$mimat}}+=$array[3];
		$hp=$array[0];
		$hp=~s/$species-//;
		$hairpin_name{$human_readable_name{$mimat}}=$hp;
		push(@ified_miRNAs,$human_readable_name{$mimat});
	    }
	}
    }
    close(DATA);
}

%seen=();
@unified_miRNAs = do { my %seen; grep { !$seen{$_}++ } @ified_miRNAs };

#for $type (@types)
#{
#    $descr=$type;
#    $descr=~s/ /_/g;
#    open(OUT,">sorted_miRNA_expression_in_primary_tumor_vs_NAT_".$descr.".tsv")#;
#    print OUT "miRNA\tConfidence\tCase\tSample ID\tSample description\tCase description\tmiRNA rpm\n";
#    close(OUT);
#}

open(OUT,">sorted_miRNA_expression_in_primary_tumor_vs_NAT.tsv");
print OUT "miRNA\tConfidence\tCase\tSample ID\tSample description\tCase description\tmiRNA rpm\n";



#for $cancer_type (@types)
#{
#    $descr=$cancer_type." ($cancer_localization{$case{$sample}})";
#    $descr=~s/ /_/g;
#    open(OUT,">sorted_miRNA_expression_in_primary_tumor_vs_NAT_".$descr.".tsv");
#    print OUT "miRNA\tConfidence\tCase\tSample ID\tSample description\tCase description\tmiRNA rpm\n";
#    close(OUT);
#}
    
for $file_path (keys %rpm)
{
    $descr=$cancer_type{$case{$samples{$file_path}[0]}}." ($cancer_localization{$case{$samples{$file_path}[0]}})";
    $descr=~s/ /_/g;

    if (($cancer_type{$case{$samples{$file_path}[0]}} ne '') && ($cancer_localization{$case{$samples{$file_path}[0]}} ne '')) #exclude cases whose cancer type or localization is not specified
    {
#	open(OUT,">>sorted_miRNA_expression_in_primary_tumor_vs_NAT_".$descr.".tsv");
	
	for $miRNA (@unified_miRNAs)
	{
	    if (grep(/^$hairpin_name{$miRNA}$/,@high_conf_hairpins))
	    {
		$annot='high_confidence';
	    }
	    else
	    {
		$annot='lower_confidence';
	    }
	    for $sample (@{$samples{$file_path}})
	    {
		if ($rpm{$file_path}{$miRNA})
		{
		    $display=$rpm{$file_path}{$miRNA};
		}
		else
		{
		    $display=0;
		}
		print OUT "$miRNA\t$annot\t$case{$sample}\t$sample\t$sample_description{$sample}\t$cancer_type{$case{$sample}} ($cancer_localization{$case{$sample}})\t$display\n";
	    }
	}
    }
}
close(OUT);

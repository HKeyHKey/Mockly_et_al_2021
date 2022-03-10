#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter script inputs (e.g., ./Module_plots_copy_number_per_cancer_type.pl gdc_sample_sheet.2021-03-04.tsv clinical.tsv).\n";
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
	@array2=split(', ',$array[5]);
	%seen=();my @unique = do { my %seen; grep { !$seen{$_}++ } @array2 };
	if (push(@unique)>1)
	{
	    print "Problem! Several cases appear to be associated with file $file_path (please correct).\n";
	    die;
	}
	$case{$file_path}=$unique[0];
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
	### shows that every case in 'clinical.tsv' (as of Mar. 4, 2021) except two (BLGSP-71-06-00252 and BLGSP-71-06-00253) has a single (cancer type, cancer localization), and except for HCM-BROD-0199-C71, HCM-BROD-0223-C43, HCM-BROD-0227-C43, HCM-BROD-0332-C43, HCM-BROD-0339-C43, HCM-CSHL-0085-C24, HCM-CSHL-0089-C25, HCM-CSHL-0434-C24 and HTMCP-03-06-02411 (for which the doublet is due to an incomplete annotation: "Tumor, NOS", "Not Reported" in some of the entries).
	### For the two cases (BLGSP-71-06-00252 and BLGSP-71-06-00253) with potentially several cancers: looks like an annotation problem (a very poorly described tumor for "Bones of skull and face and associated joints", and a precisely-described tumor for "Hematopoietic system, NOS"). We will count them only in their precisely-described cancer:
	if (($array[1] eq 'BLGSP-71-06-00252') || ($array[1] eq 'BLGSP-71-06-00253'))
	{
	    if ($array[107] ne 'Tumor, NOS')
	    {
		$cancer_type{$array[1]}=$array[107];
		$cancer_localization{$array[1]}=$array[119];
		push(@{$cases_for_cancer_type{$array[107].' ('.$array[119].')'}},$array[1]);
	    }
	}
	else
	{
	    if (($cancer_type{$array[1]} eq '') || ($array[107] ne 'Not Reported') || ($array[107] ne 'Tumor, NOS'))
	    {
		$cancer_type{$array[1]}=$array[107];
		$cancer_localization{$array[1]}=$array[119];
		push(@{$cases_for_cancer_type{$array[107].' ('.$array[119].')'}},$array[1]);
	    }
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
open(CONF_HAIRPINS,'high_conf_hairpin.fa');
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

@list_files=`ls */*.tsv`;

foreach $f (@list_files)
{
    chomp $f;
    push(@files,$f);
    open(DATA,$f);
    $type_localization=$cancer_type{$case{$f}}." (".$cancer_localization{$case{$f}}.")";
    while(<DATA>)
    {
	chomp;
	@array=split('\t',$_);
	if (($array[1]=~/^MIR/) && ($array[1]!~/HG$/))
	{
	    if ($array[5] eq '')
	    {
		$data{$array[1]}{$case{$f}}='ND';
		++$nd{$type_localization}{$array[1]};
	    }
	    else
	    {
		if ($array[5]>2)
		{
		    $data{$array[1]}{$case{$f}}='gain';
		    ++$ga{$type_localization}{$array[1]};
		}
		else
		{
		    if ($array[5]<2)
		    {
			$data{$array[1]}{$case{$f}}='loss';
			++$lo{$type_localization}{$array[1]};
		    }
		    else
		    {
			$data{$array[1]}{$case{$f}}='diploid';
			++$di{$type_localization}{$array[1]};
		    }
		}
	    }
	}
    }
    close(DATA);
}


open(OUT,">Observed_cancer_copy_number_variation_in_".$species."_miRNAs.tsv");
print OUT "# For each cancer type+localization and for each miRNA gene: number of cases with (respectively) non-determined, gain, loss and diploid miRNA gene number are given, separated by slashes (i.e.: ND/gain/loss/diploid numbers)\n";
print OUT "miRNA\tConfidence";
for $cancer_type (@types)
{
    print OUT "\t$cancer_type";
}
print OUT "\n";
    
for $weird_name (keys %data)
{
    $hairpin_name=$weird_name;
    $hairpin_name=~s/^MIRLET/let-/;
    $hairpin_name=~s/^MIR/mir-/;
    ($hairpin_prefix,$hairpin_suffix)=($hairpin_name,$hairpin_name);
    $hairpin_prefix=~s/([a-z][a-z][a-zA-Z]-\d+).*/\1/;
    $hairpin_suffix=~s/$hairpin_prefix//;
    $hairpin_suffix=~tr/QWERTYUIOPASDFGHJKLZXCVBNM/qwertyuiopasdfghjklzxcvbnm/;
    $hairpin_suffix=~s/^([a-z]*)(\d+)/\1-\2/;
    if ($hairpin_suffix ne '')
    {
	$hairpin_name=$hairpin_prefix.$hairpin_suffix;
    }
    else
    {
	$hairpin_name=$hairpin_prefix;
    }
    if (grep(/^$hairpin_name$/,@high_conf_hairpins))
    {
	$annot='high_confidence';
    }
    else
    {
	$annot='lower_confidence';
    }
    
    print OUT "$hairpin_name\t$annot";
    for $cancer_type (@types)
    {
	print OUT "\t$nd{$cancer_type}{$weird_name}/$ga{$cancer_type}{$weird_name}/$lo{$cancer_type}{$weird_name}/$di{$cancer_type}{$weird_name}";
    }
    print OUT "\n";
}
close(OUT);

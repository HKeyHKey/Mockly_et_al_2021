#!/usr/bin/perl


if ($ARGV[2] eq '')
{
    print "Please enter script argument (TSV file from GDC, describing clinical data for the downloaded MAF files; then Sample sheet from GDC, also in TSV format; then minimal number of cases for a cancer type+localization to be analyzed; e.g., ./Module_GDC_cancer_simple_nucleotide_variation_in_miRNAs.pl clinical.tsv gdc_sample_sheet.2021-02-24.tsv 100).\n";
}
else
{
    $species='hsa';
    open(CONF_HAIRPINS,'hairpin_high_conf.fa');
    while(<CONF_HAIRPINS>)
    {
	chomp;
	if (/^>$species-/)
	{
	    s/^> *//;
	    s/ .*//;
	    push(@high_conf_hairpins,$_);
	}
    }
    close(CONF_HAIRPINS);

    open(CONF_MATURE,'mature_high_conf.fa');
    while(<CONF_MATURE>)
    {
	chomp;
	if (/^>$species-/)
	{
	    s/^> *//;
	    s/ .*//;
	    push(@high_conf_mature,$_);
	}
    }
    close(CONF_MATURE);

    open(COORD,$species.'.gff3');
    while (<COORD>)
    {
	chomp;
	if ($_ !~ /^#/)
	{
	    @array=split('\t',$_);
	    ($chr,$type,$start,$end,$strand,$ID)=($array[0],$array[2],$array[3],$array[4],$array[6],$array[8]);
	    $name=$ID;
	    $ID=~s/^ID=//;
	    $ID=~s/;.*//;
	    $name=~s/.*;Name=//;
	    $name=~s/;.*//;
	    $length{$name}=$end-$start+1;
	    if ($type eq 'miRNA')
	    {
		if ($strand eq '+')
		{
		    @seed_coord=($start+1,$start+6);
		}
		else
		{
		    @seed_coord=($end-6,$end-1);
		}
		@{$mature{$name}}=($chr,$start,$end);
		@{$seed{$name}}=@seed_coord;
		if (grep(/^$name$/,@high_conf_mature))
		{
		    $annot{$name}='high_confidence';
		}
		else
		{
		    $annot{$name}='lower_confidence';
		}
	    }
	    else
	    {
		@{$hairpin{$name}}=($chr,$start,$end);
		if (grep(/^$name$/,@high_conf_hairpins))
		{
		    $annot{$name}='high_confidence';
		}
		else
		{
		    $annot{$name}='lower_confidence';
		}
	    }   
	}
    }
    close(COORD);

    print "Now parsing clinical data file...\n";

    open(CLINICAL,$ARGV[0]);
    while(<CLINICAL>)
    {
	chomp;
	if ($_ !~ /^case_id	case_submitter_id	project_id	age_at_index/)
	{
	    ### The following command:
	    # for case in `tail -n +2 clinical.tsv | awk -F '\t' '{print $2,$108,$120}' | sort | uniq | awk '{print $1}' | sort | uniq -c | grep -v '^ *1 ' | awk '{print $2}'`;do tail -n +2 clinical.tsv | awk -F '\t' '{print $2,$108,$120}' | grep -w $case;echo "";done
	    ### shows that every case in 'clinical.tsv' (as of Feb. 24, 2021) has a single (cancer type, cancer localization), except for 8 cases (HCM-BROD-0199-C71, HCM-BROD-0223-C43, HCM-BROD-0227-C43, HCM-BROD-0332-C43, HCM-BROD-0339-C43, HCM-CSHL-0085-C24, HCM-CSHL-0089-C25, HCM-CSHL-0434-C24), for which the doublet is due to an incomplete annotation (every entry except 1 have a detailed cancer type, the last one is annotated as "Not Reported")
	    ### so we can simply merge them all, and use the most precise annotation for each (i.e.: ignore the "Not Reported" entries):
	    @array=split('\t',$_);
	    if (($cancer_type{$array[1]} eq '') || ($array[107] ne 'Not Reported'))
	    {
		$cancer_type{$array[1]}=$array[107];
		$cancer_localization{$array[1]}=$array[119];
		push(@{$cases_for_cancer_type{$array[107]."\t".$array[119]}},$array[1]);
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
	
	
#    for $case (keys %cancer_type)
#    {
#	print "$case $cancer_type{$case} $cancer_localization{$case}\n";
#    }

    print "Now parsing sample sheet...\n";
    
    open(SHEET,$ARGV[1]);
    while(<SHEET>)
    {
	chomp;
	if ($_ !~ /^File ID\tFile Name\tData Category\tData Type\tProject ID\tCase ID\tSample ID\tSample Type$/)
	{
	    @array=split('\t',$_);
	    $data_file=$array[0].'/'.$array[1];
	    $data_file=~s/\.gz$//;
	    $cases=$array[5];
	    @case_array=split(', ',$cases);
	    %seen=();
	    my @unique_cases = do { my %seen; grep { !$seen{$_}++ } @case_array };
	    @{$cases_in_{$data_file}}=@unique_cases;
	    push(@all_cases,@unique_cases);
	}
    }
    close(SHEET);

    %seen=();
    my @unified = do { my %seen; grep { !$seen{$_}++ } @all_cases };
    @all_cases=@unified;

#    print "all_cases=@all_cases\n";
    
    print "Now parsing MAF files...\n";
    
    for $data_file (keys %cases_in_)
    {
#	print "Before: data_file=$data_file cases=@{$cases_in_{$data_file}}\n";
	open(MAF,$data_file);
	while(<MAF>)
	{
	    chomp;
#	    if (($_ !~ /^#/) && ($_ !~ /^Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\t/))
	    if (/^MIR/) # selects lines describing miRNA genes
	    {
		@array=split('\t',$_);
		($gene,$chr,$start,$end,$aliquot)=($array[0],$array[4],$array[5],$array[6],$array[15]);

		if (push(@{$cases_in_{$data_file}})==1) # if there is only 1 case in that MAF file, create a mock aliquot name, simply adding an underscore at the end of that case's name
		{
		    $aliquot=${$cases_in_{$data_file}}[0].'_';
		}
		push(@{$weird_names{$aliquot}},$gene);
		push(@{$var_chrs{$aliquot}},$chr);
		push(@{$var_starts{$aliquot}},$start);
		push(@{$var_ends{$aliquot}},$end);
		push(@all_aliquots,$aliquot);
	    }
	}
	close(MAF);
#	print "data_file=$data_file cases=@{$cases_in_{$data_file}} last aliquot=$aliquot\n";
    }

    %seen=();
    my @unified = do { my %seen; grep { !$seen{$_}++ } @all_aliquots };
    @all_aliquots=@unified;

#    print "all_aliquots=@all_aliquots\n";
    
    for $case (@all_cases)
    {
	@{$aliquot_for_case{$case}}=grep(/^$case[-_]/,@all_aliquots);
#	print "for case: $case aliquots=@{$aliquot_for_case{$case}}\n";
	$aliquot=${$aliquot_for_case{$case}}[0];
	if ($aliquot ne '')
	{
	#	($weird_name,$chr,$var_start,$var_end)=split(' ',$_);
	    for ($n=0;$n<push(@{$weird_names{$aliquot}});++$n)
	    {
		($chr,$var_start,$var_end)=(${$var_chrs{$aliquot}}[$n],${$var_starts{$aliquot}}[$n],${$var_ends{$aliquot}}[$n]);
		($mature_name,$hairpin_name)=(${$weird_names{$aliquot}}[$n],${$weird_names{$aliquot}}[$n]);
		$mature_name=~s/^MIRLET/let-/;
		$mature_name=~s/^MIR/miR-/;
		($mature_prefix,$mature_suffix)=($mature_name,$mature_name);
		$mature_prefix=~s/([a-z][a-z][a-zA-Z]-\d+).*/\1/;
		$mature_suffix=~s/$mature_prefix//;
		$mature_suffix=~tr/QWERTYUIOPASDFGHJKLZXCVBNM/qwertyuiopasdfghjklzxcvbnm/;
		$mature_suffix=~s/^([a-z]*)(\d*)/\1/;
		if (($mature_suffix ne '') && ($mature_suffix !~/^-\d/))
		{
		    $mature_name=$mature_prefix.$mature_suffix;
		}
		else
		{
		    $mature_name=$mature_prefix;
		}
		
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
		
		for $miRBase_mature (keys %mature)
		{
		    $stripped=$miRBase_mature;
		    $stripped=~s/-[35]p$//;
		    $stripped=~s/$species-//;
		    if ($stripped eq $mature_name)
		    {
			if (($chr eq ${$mature{$miRBase_mature}}[0]) && (($var_start<=${$mature{$miRBase_mature}}[1] && $var_end>=${$mature{$miRBase_mature}}[1]) || ($var_start<=${$mature{$miRBase_mature}}[2] && $var_end>=${$mature{$miRBase_mature}}[2]) || ($var_start<=${$mature{$miRBase_mature}}[1] && $var_end>=${$mature{$miRBase_mature}}[2]) || ($var_start>=${$mature{$miRBase_mature}}[1] && $var_end<=${$mature{$miRBase_mature}}[2])))
			{
			    ++$hit_mature{$miRBase_mature}{$cancer_type{$case}."\t".$cancer_localization{$case}};
			}
			if (($chr eq ${$mature{$miRBase_mature}}[0]) && (($var_start<=${$seed{$miRBase_mature}}[0] && $var_end>=${$seed{$miRBase_mature}}[0]) || ($var_start<=${$seed{$miRBase_mature}}[1] && $var_end>=${$seed{$miRBase_mature}}[1]) || ($var_start<=${$seed{$miRBase_mature}}[0] && $var_end>=${$seed{$miRBase_mature}}[1]) || ($var_start>=${$seed{$miRBase_mature}}[0] && $var_end<=${$seed{$miRBase_mature}}[1])))
			{
			    ++$hit_seed{$miRBase_mature}{$cancer_type{$case}."\t".$cancer_localization{$case}};
			}
		    }
		}
		
		for $miRBase_hairpin (keys %hairpin)
		{
		    $stripped=$miRBase_hairpin;
		    $stripped=~s/$species-//;
		    if ($stripped eq $hairpin_name)
		    {
			if (($chr eq ${$hairpin{$miRBase_hairpin}}[0]) && (($var_start<=${$hairpin{$miRBase_hairpin}}[1] && $var_end>=${$hairpin{$miRBase_hairpin}}[1]) || ($var_start<=${$hairpin{$miRBase_hairpin}}[2] && $var_end>=${$hairpin{$miRBase_hairpin}}[2]) || ($var_start<=${$hairpin{$miRBase_hairpin}}[1] && $var_end>=${$hairpin{$miRBase_hairpin}}[2]) || ($var_start>=${$hairpin{$miRBase_hairpin}}[1] && $var_end<=${$hairpin{$miRBase_hairpin}}[2])))
			{
			    ++$hit_hairpin{$miRBase_hairpin}{$cancer_type{$case}."\t".$cancer_localization{$case}};
			}
		    }
		}
	    }
	}
    }



    

    open(OUT,">Observed_cancer_simple_nucleotide_variation_in_".$species."_miRNAs_cutoff_".$ARGV[2].".tsv");
    print OUT "miRNA\tRecorded variations in hairpin per nt and per cancer type and localization\tRecorded variations in mature miRNA per nt and per cancer type and localization\tRecorded variations in seed per nt and per cancer type and localization\tConfidence\tCancer type\tCancer localization\tRecorded variations in hairpin\tRecorded variations in mature miRNA\tRecorded variations in seed\tNumber of cases for that cancer type and localization\n";
    for $miRBase_mature (keys %mature)
    {
	@sources=();
	@hairpin_lengths=();
	$nb=0;
	$stripped=$miRBase_mature;
	$stripped=~s/-[35]p$//;
	$stripped=~s/-miR-/-mir-/;
	for $miRBase_hairpin (keys %hairpin)
	{
	    if (($stripped eq $miRBase_hairpin) || ($miRBase_hairpin=~/^$stripped-\d+$/))
	    {
		$nb=push(@sources,$miRBase_hairpin);
		push(@hairpin_lengths,$length{$miRBase_hairpin});
	    }
	}
	if ($nb==0)
	{
	    for $miRBase_hairpin (keys %hairpin)
	    {
		$stripped_h=$miRBase_hairpin;
		$stripped_h=~s/[a-z]$//;
		if (($stripped eq $stripped_h) || ($stripped_h=~/^$stripped-\d+$/))
		{
		    $nb=push(@sources,$miRBase_hairpin);
		    push(@hairpin_lengths,$length{$miRBase_hairpin});
		}
	    }
	}

	
	$d0=0;
	for $type_localization (@types)
	{
	    if ($number_of_cases{$type_localization} >= $ARGV[2]) # selects types and localization of cancers for which a large number of cases can be analyzed
	    {
		$d0=0;
		$e0=0;
		for $h (@sources)
		{
		    $d0+=$hit_hairpin{$h}{$type_localization}/$length{$h}/$number_of_cases{$type_localization};
		    $e0+=$hit_hairpin{$h}{$type_localization};
		}
		$d0=$d0/$nb;
		$e0=$e0/$nb;
		$d1=$hit_mature{$miRBase_mature}{$type_localization}/$length{$miRBase_mature}/$number_of_cases{$type_localization};
		$e1=$hit_mature{$miRBase_mature}{$type_localization};
		$d2=$hit_seed{$miRBase_mature}{$type_localization}/6/$number_of_cases{$type_localization};
		$e2=$hit_seed{$miRBase_mature}{$type_localization};
		if ($d0 eq '')
		{
		    $d0=0;
		    $e0=0;
		}
		if ($d1 eq '')
		{
		    $d1=0;
		}
		if ($d2 eq '')
		{
		    $d2=0;
		}
		if ($e0 eq '')
		{
		    $e0=0;
		}
		if ($e1 eq '')
		{
		    $e1=0;
		}
		if ($e2 eq '')
		{
		    $e2=0;
		}
		print OUT "$miRBase_mature\t$d0\t$d1\t$d2\t$annot{$miRBase_mature}\t$type_localization\t$e0\t$e1\t$e2\t$number_of_cases{$type_localization}\n";
	    }
	}
    }
    close(OUT);
}

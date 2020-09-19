#!/bin/env perl 

use lib $ENV{"JUMP_L_LIB"};
use Getopt::Long;
use Spiders::Search;
use Spiders::Deisotope;
use Spiders::Consolidation;
use Spiders::Dta;
use Spiders::MassUtils;
#use Algorithm::Combinatorics qw/combinations/;
	
my $search = new Spiders::Search;

my ($fraction,$scan,$peptide,$outdir);

GetOptions('fraction=s'=>\$fraction,
			'scan=s'=>\$scan,
			'peptide=s'=>\$peptide,
			'outdir=s'=>\$outdir,
			'parameter=s'=>\$parameter,
		);
		
my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();
$params->{'Mn_Mn1'} = 0.5;
$params->{'M_M1'} = 0.3;


## output file
my $outfile = $scan . ".out";
open(OUTFILE,">$outdir/$outfile") || die "can not open the file: $outfile";

my $scan_num = $scan;
my @data = split(/\./,$scan_num);

# one scan 
my $dtafile = "$outdir/$scan.dta";

my $dta = new Spiders::Dta();
$dta->set_dta_file($dtafile);
$dta->process_dtafile();
$dta->parse_dtafile_name();	
my $charge = $dta->get_charge();
my $exp_mz_orig = $dta->get_mz_array();
my %exp_mz_hash_orig;
foreach (@$exp_mz_orig)
{
	my $mz_int  = int($_);
	$exp_mz_hash_orig{$mz_int}{$_}=1;
}

=head
if($params->{'MS2_deisotope'}==1)
{
	my $deisotope = new Spiders::Deisotope();
	$deisotope->set_parameter($params);
	
#	my \$dtafile = abs_path(\$dtafile);
	
	$deisotope->set_dta($dtafile);
	($ms2_signal_noise_ratio,$pho_neutral_loss) = $deisotope->MS2_deisotope();

	undef $deisotope;
}	
my $consolidation = new Spiders::Consolidation('-dta_file'=>$dtafile,'-keepnum'=>$params->{'ms2_consolidation'});
$consolidation->set_parameter($params);
$consolidation->Consolidation();
=cut
#my $peptide = $frac_scan->{$fraction}->{$scan}->{'peptide'};
my $searched_pep_orig = $peptide;

$peptide =~s/[^a-zA-Z\.]//g;
my $mutils = new Spiders::MassUtils();
$mutils->set_parameter($params);
my $peptide_mass = $mutils->get_peptide_mass_nomod($peptide);
# generate all candidate peptides
my $modif_masshash;
push (@{$modif_masshash->{$peptide_mass}}, $peptide);

#generate_modif_peptide($modif_masshash);
#my $group = groupCombinations($searched_pep_orig);

my %pvalue_hash;
print OUTFILE $fraction,"\t",$scan,"\t",$searched_pep_orig ,"\t";

## updated dtafile
$dta->set_dta_file($dtafile);
$dta->process_dtafile();
$dta->parse_dtafile_name();	
my $charge = $dta->get_charge();
my $exp_mz = $dta->get_mz_array();

my %mod_pos;
my $temp_searched_pep_orig = $searched_pep_orig;
$temp_searched_pep_orig = (split(/\./,$temp_searched_pep_orig))[1];
$temp_searched_pep_orig =~s/[^a-zA-Z]+//g;

my @seq_array = split(//,$temp_searched_pep_orig);		
for(my $i=0;$i<=$#seq_array;$i++)
{
	if(($seq_array[$i] eq 'S') or ($seq_array[$i] eq 'T') or ($seq_array[$i] eq 'Y'))
	{
		$mod_pos{$i+3}=1;
	}
}

my $modChars = "#%*";
my $groupall = allCombinations($searched_pep_orig,$modChars);
my $number_candidate_peptide = scalar keys (%$groupall);
my $site_jscore;

if($number_candidate_peptide==1)
{
#	$temp_searched_pep_orig = $searched_pep_orig;

	my @peptide_array = split(//,$temp_searched_pep_orig);
	my $number = 0;
	$number =() = $temp_searched_pep_orig =~/[^a-zA-Z]/gi;		
	for (my $i=3;$i<=$#peptide_array+3;$i++)
	{
		if($mod_pos{$i}==1)
		{		
			print OUTFILE $peptide_array[$i-3],$i-2,"\t";
			print OUTFILE 100,"\t";				
		}
	}
	print OUTFILE "\n";			
	exit;
}
else
{
	my $pep_score_total = 0;
## calculate Jscore for each peptide	
	foreach my $key (keys %$groupall)
	{
		my $pep = $groupall->{$key};
		my $mod = generate_mod($pep);		
		$search->set_parameter($params);
		my $orig_pep = $pep;
		my @pep_array = split(/\./,$orig_pep);
		$orig_pep = $pep_array[1];
		$orig_pep =~s/[^a-zA-Z]//g;
		
		$search->set_exp_mz($exp_mz);
		$search->set_precCharge($charge);
		
		my ($theo_mz,$theo_mz_hash) = $search->generate_peptide_theoretical_mass($orig_pep,$mod);
		$search->set_frag_tolerance($params->{'frag_mass_tolerance'});
		
		$exp_mz_hash = $search->get_exp_mz_hash();
				
		my ($matched, $matched_peaks) = $search->compare_theoritical_experiment(\%exp_mz_hash_orig,$theo_mz);
		my $theo_num = $#$theo_mz + 1; 
		my $pep_pvalue = sprintf("%0.40f",$search->get_peptide_pvalue($matched, $theo_num));
		if($pep_pvalue==0)
		{
			$pep_pvalue = 1E-40;
		}
		$pvalue_hash{$pep}{'score'} = 1 / $pep_pvalue;
		$pvalue_hash{$pep}{'mod'} = $mod;
		$pvalue_hash{$pep}{'orig_pep'} = $orig_pep;

		
		$pep_score_total+=$pvalue_hash{$pep}{'score'};
		print $pep,"\t",$pep_pvalue,"\t", $pvalue_hash{$pep}{'score'},"\t", $pep_score_total,"\n";
	}
## calculate the percentage of each peptide	
	foreach my $key (keys %$groupall)
	{
		my $pep = $groupall->{$key};
		$pvalue_hash{$pep}{'percentage'} = $pvalue_hash{$pep}{'score'} / $pep_score_total * 100;		
	}
##	calculate site score
	foreach my $pos (sort {$a<=>$b} keys %mod_pos)
	{		
		foreach my $key (keys %$groupall)
		{	
			my $pep = $groupall->{$key};		
			my @number = split(/\_/,$key);
			for(my $i=0;$i<=$#number;$i++)
			{
				if($number[$i] == $pos)
				{
					$site_jscore->{$pos} += $pvalue_hash{$pep}{'percentage'};
				}	
			}
		}			
	}
	foreach my $pos (sort {$site_jscore->{$b}<=>$site_jscore->{$a}} keys %$site_jscore)
	{	
		print $seq_array[$pos-3],$pos-2,"\t",sprintf("%0.2f",$site_jscore->{$pos}),"\n";
		print OUTFILE $seq_array[$pos-3],$pos-2,"\t",sprintf("%0.2f",$site_jscore->{$pos}),"\t";		
	}
	print OUTFILE "\n";
}


=head	
my $site_jscore;
foreach my $group_id (keys %$group)
{
	if($#{$group->{$group_id}}==0)
	{

	
		my @peptide_array = split(//,$temp_searched_pep_orig);
		my $number = 0;
		$number =() = $temp_searched_pep_orig =~/[^a-zA-Z]/gi;

		for (my $i=0;$i<=$#peptide_array;$i++)
		{
			if($mod_pos{$i}==1)
			{
				print OUTFILE $searched_pep_orig,"\t";			
				print OUTFILE $peptide_array[$i-1],$i-$number,"\t";
				print OUTFILE 1000,"\t";				
			}
		}
		print OUTFILE "\n";			
		exit;
	}

	my $groupall = groupAll($searched_pep_orig);

	foreach my $pos (sort {$a<=>$b} keys %mod_pos)
	{
		foreach my $key (keys %$groupall)
		{
			my $match_key = "F" . $pos . "_";
			if($key =~/$match_key/)
			{
				my $pep = $groupall->{$key};

				my $mod = generate_mod($pep);
				
				$search->set_parameter($params);
				my $orig_pep = $pep;
				my @pep_array = split(/\./,$orig_pep);
				$orig_pep = $pep_array[1];
				$orig_pep =~s/[^a-zA-Z]//g;
				
				$search->set_exp_mz($exp_mz);
				my ($theo_mz,$theo_mz_hash) = $search->generate_peptide_theoretical_mass($orig_pep,$mod);
				$search->set_frag_tolerance($params->{'frag_mass_tolerance'});
				
				$exp_mz_hash = $search->get_exp_mz_hash();
				my ($matched, $matched_peaks) = $search->compare_theoritical_experiment($exp_mz_hash,$theo_mz);
				my $theo_num = $#$theo_mz + 1; 
				my $pep_pvalue = sprintf("%0.2f",$search->get_peptide_pvalue($matched, $theo_num));
				
				$pvalue_hash{$pos}{$pep}{'pvalue'} = $pep_pvalue;
				$pvalue_hash{$pos}{$pep}{'mod'} = $mod;
				$pvalue_hash{$pos}{$pep}{'orig_pep'} = $orig_pep;
				$site_jscore->{$pos} += $pep_pvalue;			
			}
		}
	}
	if((scalar keys %$mod_pos) ==  )
}	

	
	foreach (@{$group->{$group_id}})
	{
		my $pep = $_;

		my $mod = generate_mod($pep);
		
		$search->set_parameter($params);
		my $orig_pep = $pep;
		my @pep_array = split(/\./,$orig_pep);
		$orig_pep = $pep_array[1];
		$orig_pep =~s/[^a-zA-Z]//g;
		
		$search->set_exp_mz($exp_mz);
		my ($theo_mz,$theo_mz_hash) = $search->generate_peptide_theoretical_mass($orig_pep,$mod);
		$search->set_frag_tolerance($params->{'frag_mass_tolerance'});
		
		$exp_mz_hash = $search->get_exp_mz_hash();
		my ($matched, $matched_peaks) = $search->compare_theoritical_experiment($exp_mz_hash,$theo_mz);
		my $theo_num = $#$theo_mz + 1; 
		my $pep_pvalue = sprintf("%0.2f",$search->get_peptide_pvalue($matched, $theo_num));
		
		$pvalue_hash{$group_id}{$pep}{'pvalue'} = $pep_pvalue;
		$pvalue_hash{$group_id}{$pep}{'mod'} = $mod;
		$pvalue_hash{$group_id}{$pep}{'orig_pep'} = $orig_pep;
		$site_jscore->{$group_id} += $pep_pvalue;
	}
}


	my $searched_pep_mod;		
	my $searched_pep_seq;
	foreach my $group_id (sort {$site_jscore->{$b}<=>$site_jscore->{$a}} keys %$site_jscore)
	{
		if($searched_pep_seq=="")
		{
			$searched_pep_mod = $pvalue_hash{$group_id}{$searched_pep_orig}{'mod'};
			$searched_pep_seq = $pvalue_hash{$group_id}{$searched_pep_orig}{'orig_pep'};	
		}
		
		my $sec_candidate = (sort {$pvalue_hash{$group_id}{$b}{'pvalue'}<=>$pvalue_hash{$group_id}{$a}{'pvalue'}} keys %{$pvalue_hash{$group_id}})[0];
		my $cand_pep_mod = $pvalue_hash{$group_id}{$sec_candidate}{'mod'};
		my $cand_pep_seq = $pvalue_hash{$group_id}{$sec_candidate}{'orig_pep'};

		if($sec_candidate eq $searched_pep_orig)
		{
			$sec_candidate = (sort {$pvalue_hash{$group_id}{$b}{'pvalue'}<=>$pvalue_hash{$group_id}{$a}{'pvalue'}} keys %{$pvalue_hash{$group_id}})[1];
			$cand_pep_mod = $pvalue_hash{$group_id}{$sec_candidate}{'mod'};
			$cand_pep_seq = $pvalue_hash{$group_id}{$sec_candidate}{'orig_pep'};
		}				
		next if (!defined($searched_pep_seq) or length($searched_pep_seq)==0);
		next if (!defined($cand_pep_seq) or length($cand_pep_seq)==0);
		#my $searched_pep_ions;
		my ($searched_pep_ions,$searched_pep_ions_hash) = $search->generate_peptide_theoretical_mass($searched_pep_seq,$searched_pep_mod);
		my ($cand_pep_ions,$cand_pep_ions_hash) = $search->generate_peptide_theoretical_mass($cand_pep_seq,$cand_pep_mod);
		my $pep1_diff;
		my $pep2_diff;
		($pep1_diff,$pep2_diff) = difference(\@$searched_pep_ions,\@$cand_pep_ions);

	###### Get localization score				
		my ($matched_pep1, $matched_peaks_1) = $search->compare_theoritical_experiment(\%exp_mz_hash_orig,$pep1_diff);
		my ($matched_pep2, $matched_peaks_2) = $search->compare_theoritical_experiment(\%exp_mz_hash_orig,$pep2_diff);
	#	print $searched_pep_orig,"\t",$sec_candidate,"\t",$matched_pep1,"\t",join(" ",@$pep1_diff),"\t",$matched_pep2,"\t",join(" ",@$pep2_diff),"\t";
		my $theo_num = $#$pep1_diff + 1; 
		

		my $JUMPl_score = 0;
#	if($matched_pep1>=$matched_pep2)
#	{
		my ($pep1_diff_ion_hash,$pep2_diff_ion_hash) = difference_hash($searched_pep_ions_hash,$cand_pep_ions_hash);	
		my ($exp_number,$total_number) = Generate_exp_ions($exp_mz_orig,$pep1_diff_ion_hash);
	#	print "\t",$total_number,"\t",$exp_number,"\t",$theo_num,"\t",$matched_pep1;
		my $JUMPl_score = sprintf("%0.2f",get_peptide_pvalue($total_number,$exp_number, $theo_num,$matched_pep1));
		#print "\t",$searched_pep_orig,"\t",$JUMPl_score;
#	}
#	if($matched_pep1<$matched_pep2)
#	{
#		my ($pep1_diff_ion_hash,$pep2_diff_ion_hash) = difference_hash($searched_pep_ions_hash,$cand_pep_ions_hash);	
		($exp_number,$total_number)  = Generate_exp_ions($exp_mz_orig,$pep2_diff_ion_hash);
	#	print "\t",$total_number,"\t",$exp_number,"\t",$theo_num,"\t",$matched_pep2;
		my $JUMP2_score = sprintf("%0.2f",get_peptide_pvalue($total_number,$exp_number, $theo_num,$matched_pep2));
	#	print "\t",$sec_candidate,"\t",$JUMP2_score;	
		if($JUMPl_score>=$JUMP2_score)
		{
			my $AA_site = Mod_site($searched_pep_orig,$sec_candidate);		
			print OUTFILE $searched_pep_orig,"\t",$AA_site,"\t",$JUMPl_score,"\t";			
		}
		else
		{
			my $AA_site = Mod_site($sec_candidate,$searched_pep_orig);		
			print OUTFILE $sec_candidate,"\t",$AA_site,"\t",$JUMP2_score,"\t";
			$searched_pep_orig = $sec_candidate;
			$searched_pep_mod = $cand_pep_mod;				
		}
		
	}
=cut
print OUTFILE "\n";	



	
sub Mod_site
{
	my ($searched_peptide, $candidate_peptide) = @_;
	my %mod_search;
	my %mod_cand;
	my $new_site;
	while($searched_peptide =~ /[\#\%\*]+/)
	{
		$mod_search{$-[0]}=1;
		$searched_peptide =~ s/[\#\%\*]+//;
	}
	while($candidate_peptide =~ /[\#\%\*]+/)
	{
		$mod_cand{$-[0]}=1;
		$candidate_peptide =~ s/[\#\%\*]+//;
	}	
	foreach my $site (keys %mod_search)
	{
		next if(defined($mod_cand{$site}));
		$new_site = $site;
	}
	
	my $number = 0;
	$number =() = $searched_peptide =~/[^a-zA-Z]/gi;
	my @peptide_array = split(//,$searched_peptide);
	my $pos = $new_site-$number;
	my $AA_site = $peptide_array[$new_site-1] . $pos;
	return $AA_site;
}
	
sub Generate_exp_ions
{
	my ($exp_mz_orig,$pep_ions_hash) = @_;
	my $exp_number = 0;
	my $total_number = 0;
	my @local_mz_peaks;
	foreach my $ion (keys %{$pep_ions_hash})
	{
		my @sort_ions = sort {$a<=>$b} @{$pep_ions_hash->{$ion}};

		my $mass_tolerance = $params->{'frag_mass_tolerance'};
		
		my $tolerance_unit = 1;
		if(defined($params->{'frag_mass_tolerance_unit'}))
		{
			$tolerance_unit = $params->{'frag_mass_tolerance_unit'};
		}
		
		#my $exp_mz = $self->get_exp_mz();

		
		if($tolerance_unit == 2)
		{
	### use the average mass to convert the unit	
			$mass_tolerance = $mass_tolerance *($sort_ions[0]) / 1000000;
		}
		my $min_mass;
		my $max_mass;
		if($#sort_ions==0)
		{
			$min_mass = $sort_ions[0] - 50;
			$max_mass = $sort_ions[0] + 50;
		}
		else
		{
			$min_mass = $sort_ions[0]-$mass_tolerance;
			$max_mass = $sort_ions[$#sort_ions]+$mass_tolerance;		
	## expanding 2 Da in case there is only peak available
		}
		$total_number += int(($max_mass - $min_mass)/($mass_tolerance*2));
		
	
		foreach (@$exp_mz_orig)
		{
			if($_ <= $max_mass and $_ >= $min_mass)
			{
				push(@local_mz_peaks, $_);
				$exp_number++;
			}	
		}
	}
	return ($exp_number,$total_number);
}

sub generate_mod
{
	my $nomod_pep = shift;
	my $mod='';
	@data = split(/\./,$nomod_pep);
	
	$nomod_pep =$data[1];
	my %mod_pos;
	while($nomod_pep =~ /[^a-zA-Z]+/g)
	{
		$mod_pos{$-[0]}=1;
		$nomod_pep =~ s/[^a-zA-Z]+//;
	}				
	
		my @nomodseq = split(//,$nomod_pep);
	
		my $length=length($nomod_pep);

		for(my $i=0;$i<=$length;$i++)
		{
			if($mod_pos{$i+1})
			{
				$mod .= ":";
				$mod .= $nomodseq[$i];
			}
			else
			{
				$mod .= ":";
			}
		}	
	return $mod;
}

sub difference
{
	my ($a,$b)= @_;
	my %a_hash;
	my %b_hash;
	my %All;
	my @diff_a;
	my @diff_b;
	foreach my $key (@$a,@$b) 
	{
		$All{$key}=1; 
	}
	
	foreach my $key (@$a) 
	{
		$a_hash{$key}=1; 
	}
	foreach my $key (@$b) 
	{
		$b_hash{$key}=1; 
	}
	
	foreach $key (keys %All)
	{
		if(!defined($a_hash{$key}))
		{
			push @diff_b, $key;
		}
		if(!defined($b_hash{$key}))
		{
			push @diff_a, $key;
		}		
	}
	return (\@diff_a, \@diff_b); 
}

sub difference_hash
{
	my ($a,$b)= @_;
	# my %a_hash;
	# my %b_hash;
	# my %All;
	# my @diff_a;
	# my @diff_b;
	my %diff_a_hash;
	my %diff_b_hash;
	foreach my $ion (keys %$a)
	{
		my %a_hash;
		my %b_hash;
		my %All;
		my @diff_a;
		my @diff_b;	
		foreach my $key (@{$a->{$ion}},@{$b->{$ion}}) 
		{
			$All{$key}=1; 
		}
		
		foreach my $key (@{$a->{$ion}}) 
		{
			$a_hash{$key}=1; 
		}
		foreach my $key (@{$b->{$ion}}) 
		{
			$b_hash{$key}=1; 
		}
		
		foreach $key (keys %All)
		{
			if(!defined($a_hash{$key}))
			{
				#push @diff_b, $key;
				push(@{$diff_b_hash{$ion}},$key); 
			}
			if(!defined($b_hash{$key}))
			{
	#			push @diff_a, $key;
				push(@{$diff_a_hash{$ion}},$key); 				
			}		
		}
	#	@{$diff_a_hash->{$ion}} = @diff_a;
	#	@{$diff_b_hash->{$ion}} = @diff_b;		
	}	
	return (\%diff_a_hash, \%diff_b_hash); 
}




sub generate_modif_peptide
{
	my ($pephash) = @_;
	my $cmdatabase = new Spiders::DatabaseModif();
#	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");
	my %mod_symbol = ('M'=>"@",'S'=>"#",'T'=>"%",'Y'=>"*",'G'=>"^",'K'=>"&",'D'=>'?','A'=>'~','Q'=>'!','P'=>"(",'E'=>")",'E'=>"{",'V'=>"}","V"=>"[","H"=>"]","C"=>":","F"=>",","I"=>';',"L"=>',',"R"=>"<","N"=>">","W"=>"'");	

	#	my $params = $self->get_parameter();
	
	
#	my $combNUM = $params->{'max_modif_num'};
	my $p = new Spiders::Params();
	my ($modif,$largest_modif) =$p->get_dynamic_modifications($params);

	my %new_modif=();
	my $i=0;
	foreach $modif_aa (sort keys %$modif)
	{
#		$new_modif{$modif_aa} = $modif_aa . $mod_symbol[$i];
		$new_modif{$modif_aa} = $modif_aa . $mod_symbol{$modif_aa};
		$i++;		
	}
	$cmdatabase->GenerateModifSeq($modif,$params,$pephash,\%new_modif);
	return \%new_modif;
	
}

sub groupCombinations {
       ## Input arguments
       ## $peptide = peptide sequence with modification symbols
       ## $modChars = modification symbols (e.g. $modChars = "#%*")
       
       my ($peptide) = @_;
	   my $modChars = "#%*";	   
       my @elems = split(/\./, $peptide);
       my $previousAA = $elems[0];
       my $nextAA = $elems[-1];
       $peptide = $elems[1]; 
       
       my %refModPositions;
       my $nMods = 0;
       while ($peptide =~ /([$modChars])/g) {
              $nMods++;
              ## %refModPositions
              ## Key: modification-position of the reference peptide
              ## Value: order STYs in the reference peptide
              ## e.g. peptide = "AVNS#PVNSEHKT%QLTPAAS"
              ##            There are five STYs and the modification occurs at first and third STYs
              ##            The positions of the modified STYs are 4 and 12th AA
              ##            Then, $refModPositions{4} = 1, and
              ##                     $refModPositions{12} = 3 
              ## At this step, only keys for %refModPositions are defined
              ## values are null
              $refModPositions{(pos $peptide) - $nMods} = 0;
       }
       $peptide =~ s/[$modChars]//g;
       my $nSTYs;
       my @STYs;
       while ($peptide =~ /([STY])/g) {
              $nSTYs++;
              push (@STYs, (pos $peptide));
              if (defined $refModPositions{pos $peptide}) {
                     ## Values of %refModPositions are defined here
                     $refModPositions{pos $peptide} = $nSTYs;
              }
       }      
       
       my %group;
       my $nGroup = 0;
       ## Only one modification
       if ($nMods == 1) {
              for (my $i = 0; $i < scalar(@STYs); $i++) {
                     my $groupPeptide = $peptide;
                     my $AA = substr($groupPeptide, $STYs[$i] - 1, 1);
                     if ($AA eq "S") {
                           substr ($groupPeptide, $STYs[$i], 0, "#");
                     } elsif ($AA eq "T") {
                           substr ($groupPeptide, $STYs[$i], 0, "%");
                     } elsif ($AA eq "Y") {
                           substr ($groupPeptide, $STYs[$i], 0, "*");
                     } else {
                           die "Modification defined at the AA other than STY\n";
                     }
                     $groupPeptide = $previousAA.".".$groupPeptide.".".$nextAA;
                     push (@{$group{$nGroup}}, $groupPeptide);
              }
       ## Multiple modifications
       } elsif ($nMods == 2) {
              foreach my $key (keys %refModPositions) {
                     my @v = @STYs;
                     splice (@v, $refModPositions{$key} - 1, 1);
                     my @c = combinations(\@v, $nMods - 1);
                     for (my $i = 0; $i < scalar(@c); $i++) {
                           my $groupPeptide =$peptide;
                           $c[$i] = ([$STYs[$refModPositions{$key} - 1], @{$c[$i]}]);
                           @{$c[$i]} = sort {$a <=> $b} (@{$c[$i]});              
                           for (my $j = 0; $j < scalar(@{$c[$i]}); $j++) {
                                  my $AA = substr($groupPeptide, $c[$i][$j] + $j - 1, 1);
                                  if ($AA eq "S") {
                                         substr ($groupPeptide, $c[$i][$j] + $j, 0, "#");
                                  } elsif ($AA eq "T") {
                                         substr ($groupPeptide, $c[$i][$j] + $j, 0, "%");
                                  } elsif ($AA eq "Y") {
                                         substr ($groupPeptide, $c[$i][$j] + $j, 0, "*");
                                  } else {
                                         die "Modification defined at the AA other than STY\n";
                                  }
                           }
                           $groupPeptide = $previousAA.".".$groupPeptide.".".$nextAA;
                           push (@{$group{$nGroup}}, $groupPeptide);              
                     }
                     $nGroup++;
              }
       } elsif ($nMods == 3) {
              ## Fix one modification site
              # foreach my $key (keys %refModPositions) {
                     # my @v = @STYs;
                     # splice (@v, $refModPositions{$key} - 1, 1);
                     # my @c = combinations(\@v, $nMods - 1);
                     # for (my $i = 0; $i < scalar(@c); $i++) {
                           # my $groupPeptide =$peptide;
                           # $c[$i] = ([$STYs[$refModPositions{$key} - 1], @{$c[$i]}]);
                           # @{$c[$i]} = sort {$a <=> $b} (@{$c[$i]});              
                           # for (my $j = 0; $j < scalar(@{$c[$i]}); $j++) {
                                  # my $AA = substr($groupPeptide, $c[$i][$j] + $j - 1, 1);
                                  # if ($AA eq "S") {
                                         # substr ($groupPeptide, $c[$i][$j] + $j, 0, "#");
                                  # } elsif ($AA eq "T") {
                                         # substr ($groupPeptide, $c[$i][$j] + $j, 0, "%");
                                  # } elsif ($AA eq "Y") {
                                         # substr ($groupPeptide, $c[$i][$j] + $j, 0, "*");
                                  # } else {
                                         # die "Modification defined at the AA other than STY\n";
                                  # }
                           # }
                           # $groupPeptide = $previousAA.".".$groupPeptide.".".$nextAA;
                           # push (@{$group{$nGroup}}, $groupPeptide);              
                     # }
                     # $nGroup++;
              # }
              ## Fix two modification sites
              my @modPositions = keys %refModPositions;
              my @fixedPositions = combinations(\@modPositions, $nMods - 1);       
              foreach my $fixedPosition (@fixedPositions) {
                     @{$fixedPosition} = sort {$a <=> $b} @{$fixedPosition};
                     my @v = @STYs;
                     for (my $i = 0; $i < 2; $i++) {
                           splice (@v, $refModPositions{$$fixedPosition[$i]} - $i - 1, 1);
                     }
                     my @c;
                     for (my $i = 0; $i < scalar(@v); $i++) {
                           my $groupPeptide = $peptide;
                           $c[$i] = ([@{$fixedPosition}, $v[$i]]);
                           @{$c[$i]} = sort {$a <=> $b} (@{$c[$i]});
                           for (my $j = 0; $j < scalar(@{$c[$i]}); $j++) {
                                  my $AA = substr($groupPeptide, $c[$i][$j] + $j - 1, 1);
                                  if ($AA eq "S") {
                                         substr ($groupPeptide, $c[$i][$j] + $j, 0, "#");
                                  } elsif ($AA eq "T") {
                                         substr ($groupPeptide, $c[$i][$j] + $j, 0, "%");
                                  } elsif ($AA eq "Y") {
                                         substr ($groupPeptide, $c[$i][$j] + $j, 0, "*");
                                  } else {
                                         die "Modification defined at the AA other than STY\n";
                                  }
                           }
                           $groupPeptide = $previousAA.".".$groupPeptide.".".$nextAA;
                           push (@{$group{$nGroup}}, $groupPeptide);
                     }
                     $nGroup++;
              }
       } else {
              die "Cannot deal with more than 3 modification sites\n";
       }
       return (\%group);
}

sub groupAll {
       my ($peptide) = @_;
#      my $peptide = "K.AAS#ASAAT#AAASAAAT#AAA.R";
      my $modChars = "#%*";
       my @elems = split(/\./, $peptide);
       my $previousAA = $elems[0];
       my $nextAA = $elems[-1];
       $peptide = $elems[1]; 
       
       my %refModPositions;
       my $nMods = 0;
       while ($peptide =~ /([$modChars])/g) {
              $nMods++;
              ## %refModPositions
              ## Key: modification-position of the reference peptide
              ## Value: order STYs in the reference peptide
              ## e.g. peptide = "AVNS#PVNSEHKT%QLTPAAS"
              ##            There are five STYs and the modification occurs at first and third STYs
              ##            The positions of the modified STYs are 4 and 12th AA
              ##            Then, $refModPositions{4} = 1, and
              ##                     $refModPositions{12} = 3 
              ## At this step, only keys for %refModPositions are defined
              ## values are null
              $refModPositions{(pos $peptide) - $nMods} = 0;
       }
       $peptide =~ s/[$modChars]//g;
       my $nSTYs;
       my @STYs;
       while ($peptide =~ /([STY])/g) {
              $nSTYs++;
              push (@STYs, (pos $peptide));
              if (defined $refModPositions{pos $peptide}) {
                     ## Values of %refModPositions are defined here
                     $refModPositions{pos $peptide} = $nSTYs;
              }
       }      
       
       my %comb;
       ## Only one modification
       if ($nMods == 1) {
              for (my $i = 0; $i < scalar(@STYs); $i++) {
                     my $combPeptide = $peptide;
                     my $combKey = "V".($STYs[$i] + 2); ## +2 => count the first peptide (K/R) and a dot
                     my $AA = substr($combPeptide, $STYs[$i] - 1, 1);
                     if ($AA eq "S") {
                           substr ($combPeptide, $STYs[$i], 0, "#");
                     } elsif ($AA eq "T") {
                           substr ($combPeptide, $STYs[$i], 0, "%");
                     } elsif ($AA eq "Y") {
                           substr ($combPeptide, $STYs[$i], 0, "*");
                     } else {
                           die "Modification defined at the AA other than STY\n";
                     }
                     $combPeptide = $previousAA.".".$combPeptide.".".$nextAA;
                     $comb{$combKey} = $combPeptide;
              }
       ## Multiple modifications
       } elsif ($nMods == 2) {
              foreach my $key (keys %refModPositions) {
                     my @v = @STYs;
                     splice (@v, $refModPositions{$key} - 1, 1);
                     my @c;
                     for (my $i = 0; $i < scalar(@v); $i++) {
                           my $combPeptide = $peptide;
                           my $combKey = "F".($key + 2)."_"."V".($v[$i] + 2);     ## +2 => count the first peptide (K/R) and a dot
                           $c[$i] = ([$key, $v[$i]]);
                           @{$c[$i]} = sort {$a <=> $b} (@{$c[$i]});              
                           for (my $j = 0; $j < scalar(@{$c[$i]}); $j++) {
                                  my $AA = substr($combPeptide, $c[$i][$j] + $j - 1, 1);
                                  if ($AA eq "S") {
                                         substr ($combPeptide, $c[$i][$j] + $j, 0, "#");
                                  } elsif ($AA eq "T") {
                                         substr ($combPeptide, $c[$i][$j] + $j, 0, "%");
                                  } elsif ($AA eq "Y") {
                                         substr ($combPeptide, $c[$i][$j] + $j, 0, "*");
                                  } else {
                                         die "Modification defined at the AA other than STY\n";
                                  }
                           }
                           $combPeptide = $previousAA.".".$combPeptide.".".$nextAA;                   
                           $comb{$combKey} = $combPeptide;          
                     }
              }
       } elsif ($nMods == 3) {
              ## Fix two modification sites
              my @modPositions = keys %refModPositions;
              my @fixedPositions = combinations(\@modPositions, $nMods - 1);       
              foreach my $fixedPosition (@fixedPositions) {
                     @{$fixedPosition} = sort {$a <=> $b} @{$fixedPosition};
                     my @v = @STYs;
                     my @keyArray;
                     for (my $i = 0; $i < 2; $i++) {
                           splice (@v, $refModPositions{$$fixedPosition[$i]} - $i - 1, 1);
                           push (@keyArray, "F".($$fixedPosition[$i] + 2));       ## +2 => count the first peptide (K/R) and a dot
                     }
                     my $combKey = join("_", @keyArray);
                     my @c;
                     for (my $i = 0; $i < scalar(@v); $i++) {
                           my $combPeptide = $peptide;
                           my $combKey = $combKey."_"."V".($v[$i] + 2);
                           $c[$i] = ([@{$fixedPosition}, $v[$i]]);
                           @{$c[$i]} = sort {$a <=> $b} (@{$c[$i]});
                           for (my $j = 0; $j < scalar(@{$c[$i]}); $j++) {
                                  my $AA = substr($combPeptide, $c[$i][$j] + $j - 1, 1);
                                  if ($AA eq "S") {
                                         substr ($combPeptide, $c[$i][$j] + $j, 0, "#");
                                  } elsif ($AA eq "T") {
                                         substr ($combPeptide, $c[$i][$j] + $j, 0, "%");
                                  } elsif ($AA eq "Y") {
                                         substr ($combPeptide, $c[$i][$j] + $j, 0, "*");
} else {
                                         die "Modification defined at the AA other than STY\n";
                                  }
                           }
                           $combPeptide = $previousAA.".".$combPeptide.".".$nextAA;
                           $comb{$combKey} = $combPeptide;
                     }
              }
       } else {
              die "Cannot deal with more than 3 modification sites\n";
       }
       return (\%comb);
}

sub allCombinations {
       my ($peptide, $modChars) = @_;
       my @elems = split(/\./, $peptide);
       my $previousAA = $elems[0];
       my $nextAA = $elems[-1];
       $peptide = $elems[1]; 
       
       ## Count the number of modifications
       my $nMods = 0;
       while ($peptide =~ /([$modChars])/g) {
              $nMods++;
       }
       
       ## Identify the position(s) of STY(s)
       $peptide =~ s/[$modChars]//g;
	   $peptide =~ s/\@//g;
	   
       my @STYs;
       while ($peptide =~ /([STY])/g) {
              push (@STYs, (pos $peptide));
       }
       
       ## Generate all possible combinations of modification and store it to %comb
       my %comb;
       my @c = combinations(\@STYs, $nMods);
       for (my $i = 0; $i < scalar(@c); $i++) {
              my $combPep = $peptide;
              for (my $j = 0; $j < scalar(@{$c[$i]}); $j++) {
                     my $AA = substr($combPep, $c[$i][$j] + $j - 1, 1);
                     if ($AA eq "S") {
                           substr ($combPep, $c[$i][$j] + $j, 0, "#");
                     } elsif ($AA eq "T") {
                           substr ($combPep, $c[$i][$j] + $j, 0, "%");
                     } elsif ($AA eq "Y") {
                           substr ($combPep, $c[$i][$j] + $j, 0, "*");
                     } else {
                           die "Modification defined at the AA other than STY\n";
                     }             
              }
              my @keyArray;
              foreach (@{$c[$i]}) {
                     push (@keyArray, $_ + 2);
              }
              my $key = join("_", @keyArray);
              $comb{$key} = $previousAA . "." . $combPep . "." . $nextAA;           
       }
       return (\%comb);
}

sub combinations {
	my ($list, $n) = @_;
	die "Insufficient list members" if $n > @$list;
	return map [$_], @$list if $n <= 1;
	my @comb;
	for (my $i = 0; $i+$n <= @$list; ++$i) {
		my $val  = $list->[$i];
		my @rest = @$list[$i+1..$#$list];
		push (@comb, [$val, @$_]) for combinations(\@rest, $n - 1);
	}
	return @comb;
}

sub get_peptide_pvalue
{
	my ($total,$exp_num,$theo_num,$matched)=@_;
######### use different method for matching score #############
	
	my $hyper = new Spiders::Hypergeometric();

	my $log_peptide_pvalue = 1;

	my $peptide_pvalue=$hyper->Hypergeometric($total,$exp_num,$theo_num,$matched);	
	$log_peptide_pvalue = sprintf("%.6f",-log($peptide_pvalue)/log(10));
	return ($log_peptide_pvalue);
}

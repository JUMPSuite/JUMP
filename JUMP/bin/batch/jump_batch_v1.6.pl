#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use Cwd;
my $cwd_dir = getcwd;

if (scalar(@ARGV)!=1)
{
        die "Usage: jump_batch jump_batch.params\n";
}

# R code path
my (%parahash);
$parahash{'R_code_path'}="$Bin/R/batch_normalization.R";
#$parahash{'R_code_path'}='/home/yli4/development/JUMPbatch/create/batch_normalization.R';

# initialization
my (%uniqID,%features,%expression,%batches,%ft2name,%featureIndex,%normexp);
parse_params($ARGV[0],\%parahash);
# initialize features
featureInitialization(\%ft2name);

# mkdir output folder
my $outputDir="batch_$parahash{output_folder}";
if (-e $outputDir) {} else { system("mkdir $outputDir"); }
system("cp $ARGV[0] $outputDir");
#my $logfile="$outputDir/jump_i.log";
#open(LOG,">$logfile");

# read Data matrix
# parse id_uni_prot_quan.txt: build unique row names
# %uniqID{$id}
print "Parsing id_uni_prot_quan.txt\n";
foreach my $batch (keys %{$parahash{'batches'}}) {

	my ($input,$idCol);
	if ($parahash{'input_mode'} == 1) {
		$input="$parahash{'batches'}{$batch}{path}/id_uni_prot_quan.txt";
		$idCol=1;
	} elsif ($parahash{'input_mode'} == 2) {
		$input="$parahash{'batches'}{$batch}{path}/id_uni_site_quan.txt";
		$idCol=1;
	} elsif ($parahash{'input_mode'} == 3) {
		$input="$parahash{'batches'}{$batch}{path}/id_uni_pep_quan.txt";
		$idCol=0;
	} elsif ($parahash{'input_mode'} == 4) {
		$input="$parahash{'batches'}{$batch}{path}/id_all_prot_quan.txt";
		$idCol=1;
	} else {
		die "unexpected value for input_mode ($parahash{'input_mode'})!!!\n";
	}

	open(IN,$input) || die "cannot open $input";
	my $line=<IN>;
	$line=<IN>;

	while (<IN>) {
		chomp;
		chop if (m/\r$/); # windows ending
		my $id=(split /\t/,$_)[$idCol];
		$id=validID($id); # convert '#' to '*'
		#next unless ($id =~ /$parahash{focus_species}/);
		$uniqID{$id}='';
	}
	close IN;
}

# parse id_all_prot_quan.txt
# %features{$id}{common}{$feature}=value
#               {batchSpecific}{$batch}{$feature}=value
# %expression{$id}{$batch}{$sample}=value
# %batches{$batch}[$k]=$sample
# %featureIndex{common}{$feature}=$col_i
#              {batchSpecific}{$feature}=$col_i
print "Parsing id_all_prot_quan.txt\n";
my $minTMTint; # to replace 'NA'
foreach my $batch (keys %{$parahash{'batches'}}) {
	my ($input,$idCol);

	if (defined($parahash{'isoform_rescue'}) and $parahash{'isoform_rescue'} == 0) { # skipp isoform rescue
		if ($parahash{'input_mode'} == 1) {
			$input="$parahash{'batches'}{$batch}{path}/id_uni_prot_quan.txt";
			$idCol=1;
		} elsif ($parahash{'input_mode'} == 2) {
			$input="$parahash{'batches'}{$batch}{path}/id_uni_site_quan.txt";
			$idCol=1;
		} elsif ($parahash{'input_mode'} == 3) {
			$input="$parahash{'batches'}{$batch}{path}/id_uni_pep_quan.txt";
			$idCol=0;
		} elsif ($parahash{'input_mode'} == 4) {
			$input="$parahash{'batches'}{$batch}{path}/id_all_prot_quan.txt";
			$idCol=1;
		} else {
			die "unexpected value for input_mode ($parahash{'input_mode'})!!!\n";
		}
	} else { 	# turn on isoform rescue
		if ($parahash{'input_mode'} == 1) {
			$input="$parahash{'batches'}{$batch}{path}/id_all_prot_quan.txt";
			$idCol=1;
		} elsif ($parahash{'input_mode'} == 2) {
			$input="$parahash{'batches'}{$batch}{path}/id_all_site_quan.txt";
			$idCol=1;
		} elsif ($parahash{'input_mode'} == 3) {
			$input="$parahash{'batches'}{$batch}{path}/id_all_pep_quan.txt";
			$idCol=0;
		} elsif ($parahash{'input_mode'} == 4) {
			$input="$parahash{'batches'}{$batch}{path}/id_all_prot_quan.txt";
			$idCol=1;
		} else {
			die "unexpected value for input_mode ($parahash{'input_mode'})!!!\n";
		}
	}

	$parahash{'input_mode'} =1 if ($parahash{'input_mode'} == 4); # afterwards, mode '4' will be treated as '1', as they are both on 'protein' level

	open(IN,$input) || die "cannot open $input";
	my $line=<IN>;
        $line=<IN>;
	parseHeader($batch,$line,\%parahash,\%batches,\%ft2name,\%featureIndex);

        while (<IN>) {
                chomp;
		chop if (m/\r$/); # windows ending
		my @t=split /\t/,$_;
		my $id=$t[$idCol];
		$id=validID($id); # convert '#' to '*'

		next unless (defined($uniqID{$id}));

		# expression
		my $n=$parahash{'batches'}{$batch}{n};
		for (my $i=0;$i<$n;$i++) {
			$expression{$id}{$batch}{$batches{$batch}[$i]}=$t[$#t-$n+1+$i];
		}

=head
deactivated on 10/11/18
reason: this is based on the assumption that every channel is quantified if PSM>0; but now this assumption is not required anymore
		# min int
		if (!defined($minTMTint) || $minTMTint>min(@t[$#t-$n+1 .. $#t])) {
			$minTMTint=min(@t[$#t-$n+1 .. $#t]);
		}
=cut
		# features
		# common
		if (!defined($features{$id}{common}{GN})) {			
			foreach my $f (keys %{$featureIndex{common}}){			
				unless (defined($features{$id}{common}{$f})) {
					$features{$id}{common}{$f}='';
				}
				if (defined($t[$featureIndex{common}{$f}])) {
					$features{$id}{common}{$f}=$t[$featureIndex{common}{$f}];
				} 
			}
		}
		# batch specific
		foreach my $f (keys %{$featureIndex{batchSpecific}}){
			if (defined($t[$featureIndex{batchSpecific}{$f}])) {
				$features{$id}{batchSpecific}{$batch}{$f}=$t[$featureIndex{batchSpecific}{$f}];
			} else {
				$features{$id}{batchSpecific}{$batch}{$f}='NA';
			}
		}
        }
        close IN;

}
#print "Min TMT intensity: $minTMTint\n";

# for proteins with n/a, calculate the mean(log10 TMT) across all batches that are quantified; then fill in the n/a with that mean
my %meanLog10TMT;
my %Om; # observation matrix: $Om{$id}{$batch}{$batches{$batch}[$i]} = 1/0;

foreach my $id (keys %expression) {
	my @tmp; # array to hold the meanLog10TMT for each batch;

	# for protein, record in how many samples, quantificaiton is missed
	$features{$id}{common}{missingN}=0; 

	foreach my $batch (keys %batches) {

		# if that batch is quantified (assuming as long as the 1st sample is quantified, the whole batch is quantified; this assumption is generally correct)
		#if (defined($expression{$id}{$batch}{$batches{$batch}[0]}) 
		#and $expression{$id}{$batch}{$batches{$batch}[0]}>0) { 
		#	push @tmp, mean(log10(@{$expression{$id}{$batch}{$batches{$batch}}}));
		#}

		my $n=$parahash{'batches'}{$batch}{n};
		my @tmp2; # array to hold log10TMT for individual sample
		for (my $i=0;$i<$n;$i++) {
			if (defined($expression{$id}{$batch}{$batches{$batch}[$i]})
			and $expression{$id}{$batch}{$batches{$batch}[$i]} ne 'NA') {
				push @tmp2, log10($expression{$id}{$batch}{$batches{$batch}[$i]});
				$Om{$id}{$batch}{$batches{$batch}[$i]} = 1; # observed
			} else {
				$Om{$id}{$batch}{$batches{$batch}[$i]} = 0; # missed
				$features{$id}{common}{missingN}++;
			}
		}

		if (scalar(@tmp2) > 0) { # if at least one sample is quantified for this batch
			push @tmp, mean(@tmp2);
		}

	}

	$meanLog10TMT{$id}=mean(@tmp);
}

# printing
if (!(-e "$outputDir/intermediate")) {
	system(qq(mkdir $outputDir/intermediate));
}
# expression table
open(EXP,">$outputDir/intermediate/raw_exp.txt");
open(BCH,">$outputDir/intermediate/batch_vector.txt");
# header
foreach my $batch (keys %batches) {
	my $n=$parahash{'batches'}{$batch}{n};
	for (my $i=0;$i<$n;$i++) {
		print EXP "\t$batch\.$batches{$batch}[$i]";
		print BCH "$batch\n";
	}
}
print EXP "\n";
# data
foreach my $id (keys %expression) {
	print EXP "$id";
	foreach my $batch (keys %batches) {
		my $n=$parahash{'batches'}{$batch}{n};
		for (my $i=0;$i<$n;$i++) {
			if (defined($expression{$id}{$batch}{$batches{$batch}[$i]}) 
			and $expression{$id}{$batch}{$batches{$batch}[$i]} ne 'NA') {
				print EXP "\t$expression{$id}{$batch}{$batches{$batch}[$i]}";
			} else {
				#if (defined($minTMTint) && $minTMTint>0) {
				#	print EXP "\t",$minTMTint/2;
				#} else {
				#	print EXP "\t0";
				#}
				if (defined($meanLog10TMT{$id}) && $meanLog10TMT{$id}>0) {
					print EXP "\t",10**$meanLog10TMT{$id};
				} else {
					print "WARNING: mean TMT intensity cannot be obtained for $id; this entry will be ignored\n";
				}
			}
		}
	}
	print EXP "\n";
}
close EXP;
close BCH;

# internal_standard_vector.txt
open(STD,">$outputDir/intermediate/internal_standard_vector.txt");
foreach my $batch (keys %batches) {
	if (defined($parahash{'batches'}{$batch}{internal_standarad}) and 
	$parahash{'batches'}{$batch}{internal_standarad} =~ /^sig/) {
#	if ($parahash{'batches'}{$batch}{internal_standarad}>0) { # when using # (old code)
		my $type=($parahash{'normalization_standard'} eq $batch)?'standard':'non_standard';
#		print STD "$batch\t$batch\.$batches{$batch}[$parahash{'batches'}{$batch}{internal_standarad}-1]\t$type\n";  # when using # (old code)
		print STD "$batch\t$batch\.$parahash{'batches'}{$batch}{internal_standarad}\t$type\n";
	} else {
		print STD "NA\n";
	}
}
close STD;

# print final publication table
print "Printing publication tables\n";
if (!(-e "$outputDir/publication")) {
	system(qq(mkdir $outputDir/publication));
}
printPublicationTable("$outputDir/publication/combined_raw_uni_prot.txt",\%expression,\%features,\%batches,$parahash{input_mode});
printPublicationTable("$outputDir/publication/combined_raw_uni_prot_R.txt",\%expression,\%features,\%batches,$parahash{input_mode},1);

# For using jump -i: 
# 1. print output table in -i input format
# 2. print jump_i.params
#print "Printing tables for using jump -i\n";
if (!(-e "$outputDir/publication/jump_i")) {
	system(qq(mkdir $outputDir/publication/jump_i));
}
printPublicationTable("$outputDir/publication/jump_i/overlapped_raw_uni_prot.txt",\%expression,\%features,\%batches,$parahash{input_mode},2);

# batch normalization
if ( $parahash{'normalization_method'} == 0 ) {
	print "Skip expression normlization\n";
} else {
	if ($parahash{'normalization_method'} == 1) {
		$parahash{'normalization_method'} = 'InternalStandard';
		print "Normalize expression table across batches using internal standards\n";
	} elsif ($parahash{'normalization_method'} == 2) {
		$parahash{'normalization_method'} = 'LinearModel';
		print "Normalize expression table across batches using linear model fitting\n";
	}

	system(qq(R CMD BATCH --no-save --args -method=$parahash{'normalization_method'} -raw_exp=$outputDir/intermediate/raw_exp.txt -batch_vector=$outputDir/intermediate/batch_vector.txt -output=$outputDir/intermediate/norm_exp.txt -internal_standard=$outputDir/intermediate/internal_standard_vector.txt $parahash{'R_code_path'} $outputDir/intermediate/R.log));

	# convert norm_exp.txt 2 hash
	norm_exp2hash("$outputDir/intermediate/norm_exp.txt",\%normexp,\%batches,\%parahash);
	# print final publication table
	print "Printing publication tables\n";
	printPublicationTable("$outputDir/publication/combined_norm_uni_prot.txt",\%normexp,\%features,\%batches,$parahash{input_mode});
	printPublicationTable("$outputDir/publication/combined_norm_uni_prot_R.txt",\%normexp,\%features,\%batches,$parahash{input_mode},1);
	printPublicationTable("$outputDir/publication/jump_i/overlapped_norm_uni_prot.txt",\%normexp,\%features,\%batches,$parahash{input_mode},2);

	print_jump_i_params("$outputDir/publication/jump_i/jump_i_custom.params","$cwd_dir/$outputDir/publication/jump_i/overlapped_norm_uni_prot.txt",\%batches);

	# run jump -i
	print "\nTrying to run -i to generate heatmaps for normalized data:\n\n";
	if (defined($parahash{jump_i_path}) and $parahash{jump_i_path} ne '0') {
		# use program specified by parameter
		system(qq(cd $outputDir/publication/jump_i; perl $parahash{jump_i_path} jump_i_custom.params));
		#my $cmd="cd $outputDir/publication/jump_i; perl $parahash{jump_i_path} jump_i_custom.params";
		#print "$cmd\n";
		#system($cmd);
	} else {
		# use program in production
		system(qq(cd $outputDir/publication/jump_i; jump -i jump_i_custom.params));
	}
}


print "\njump -q-batch is finished\n";
#--------------------------------------------------------------------------------
sub print_jump_i_params {
	my ($outputFile, $jump_i_input, $batches)=@_;

	open(OUT,">$outputFile");

	# input file path
	print OUT "\n\# customized input\n";
	print OUT "input_table = $jump_i_input\n";

	# column names
	print OUT "\n\# sample labels\n";
	#foreach my $batch (sort {$a <=> $b} keys %{$batches}) {
	my @tmp; foreach (keys %{$batches}) { s/batch//; push @tmp, $_; }
	foreach (sort {$a <=> $b} @tmp) {
		my $batch = "batch$_";
		for (my $i=0; $i<scalar(@{$$batches{$batch}});$i++) {
			print OUT "$batch\_$$batches{$batch}[$i] = $batch\_sample",$i+1,"\n";
		}
	}

	# other parameters
	print OUT "\n\# output folder\n";
	print OUT "output_folder = test1\n";
	print OUT "\n\# contaminant removel option\n";
	print OUT "remove_contaminants = 0\n";
	print OUT "\n\# SD / mean matrix: analysis based on the stable proteins / peptides\n";
	print OUT "pair_cutoff_percentage = 10\n";
	print OUT "pair_cutoff_intensity = 0\n";
	print OUT "\n\# heatmap / sample clustering: analysis based on the most variable proteins / peptides\n";
	print OUT "bypass_row_clustering = 0\n";
	print OUT "cluster_cutoff_percentage = 100\n";
	print OUT "cluster_cutoff_intensity = 0\n";

	close OUT;
}

sub validID {
	my ($id)=@_;

	$id =~ s/\#/\*/g;
	return $id;
}

sub norm_exp2hash {
	my ($input,$expression,$batches,$parahash)=@_;
	open (IN,$input) || die "cannot open $input\n";
	my $line=<IN>;
	while (<IN>) {
		chomp;
		chop if (m/\r$/); # windows ending
		my @t=split /\t/,$_;
		my $id=shift @t;
		foreach my $batch (keys %{$batches}) {
			my $n=$$parahash{'batches'}{$batch}{n};
			for (my $i=0;$i<$n;$i++) {
				$$expression{$id}{$batch}{$$batches{$batch}[$i]}=shift @t;
			}
		}
	}
	close IN;
}

sub printPublicationTable {
# input_mode: 1: protein; others (2/3): peptide or site
# output_mode: 0 (or not defined): publication table; 1: R format w/ all information; 2: -i input format (only ID column and abundance columns; only proteins w/o n/a)
	my ($output,$expression,$features,$batches,$input_mode,$output_mode)=@_;

	#open (OUT,">$output");
	open (OUT,">.tmp_body"); # print body to tmp file

	if (defined($output_mode) and $output_mode==1) { # R format w/ all information
		if ($input_mode==1) { # protein mode
			print OUT 'ProteinGroup	ProteinAccession	ProteinDescription	GN	GPCRs	TFs	epigenetic_regularos	kinases	oncogenes';
		} else { # peptide or site mode
			print OUT 'ID	ModSites	ProteinGroup	ProteinAccession	ProteinDescription	GN';
		}
	} elsif (defined($output_mode) and $output_mode==2) { # -i input format
		print OUT 'customID';
	} elsif (!defined($output_mode) or $output_mode==0) { # publication table format
		if ($input_mode==1) {
			print OUT 'Protein Group#	Protein Accession #	Protein Description	GN	GPCRs	TFs	epigenetic_regularos	kinases	oncogenes';
		} else {
			print OUT 'ID	ModSites	Protein Group#	Protein Accession #	Protein Description	GN';
		}
	}

	#foreach my $batch (sort {$a <=> $b} keys %{$batches}) {
	my @tmp; foreach (keys %{$batches}) { s/batch//; push @tmp, $_; }
	foreach (sort {$a <=> $b} @tmp) {
		my $batch = "batch$_";
		#unless (defined($output_mode) and $output_mode==2) {
		if (!defined($output_mode) or $output_mode==0 or $output_mode==1) {
			print OUT "\t$batch";
		}

		if (defined($output_mode) and $output_mode==1) {
			if ($input_mode==1) {
				print OUT "	$batch\_PSMs	$batch\_TotalPeptides	$batch\_UniquePeptides	$batch\_Coverage";
			} else {
				print OUT "\t$batch\_PSMs";
			}
			for (my $i=0; $i<scalar(@{$$batches{$batch}});$i++) {
				print OUT "\t$batch\_$$batches{$batch}[$i]";
			}
		} elsif (!defined($output_mode) or $output_mode==0) {
			if ($input_mode==1) {
				print OUT '	PSM#	Total Peptide#	Unique Peptide#	%Coverage';
			} else {
				print OUT '	PSM#';
			}
			for (my $i=0; $i<scalar(@{$$batches{$batch}});$i++) {
				print OUT "\t$$batches{$batch}[$i]";
			}
		} elsif (defined($output_mode) and $output_mode==2) {
			for (my $i=0; $i<scalar(@{$$batches{$batch}});$i++) {
				print OUT "\t$batch\_$$batches{$batch}[$i]";
			}
		}
	}
	print OUT "\n";

	my %counts; # record in how many batches this protein is ID
	my $totalCount=0;
	my %singleBatchCounts; # record batch specific protein numbers
	foreach my $batch (keys %{$batches}) {
		$singleBatchCounts{$batch}=0;
	}

	foreach my $id (keys %{$expression}) {
=head
#code deactivated due to validation of assumption (on 10/11/18)
		# check in how many batches this protein are quantified (added on 2/9/18)
		# WARNING: this calculation is based on the assumption that if for one batch, its PSMs > 0, all the TMT channels are quantified
		#my $bk0=scalar(keys %{$$features{$id}{batchSpecific}}); # wrong statement
		my $bk0=0;
		foreach my $batch (sort {$a <=> $b} keys %{$batches}) {
			if (defined($$features{$id}{batchSpecific}{$batch}{PSMs}) and
			$$features{$id}{batchSpecific}{$batch}{PSMs}>0) {
				$bk0++;
			}
		}
=cut

		my $wNA=0; 	# indicator: whether n/a exist for this protein
		#$wNA=1 if ($bk0<scalar(keys %{$batches})); # old statement
				# WARNING: this is based on the assumption that if for one batch, its PSMs > 0, all the TMT channels are quantified
				# this function should be revised if the above assumption is not valid; otherwise, the -i input will be invalid!!!

		# to solve the above issue, now we define a new entry: {common}{missingN}
		$wNA = 1 if $$features{$id}{common}{missingN} > 0; # updated on 10/11/18

		# initialization
		$totalCount++;
		my $bk=0; # count in how many batches this protein is ID
		my @IDbatches; # record in which batches this protein is ID

		# common features
		if (!defined($output_mode) or $output_mode==0 or $output_mode==1) {
			if ($input_mode==1) {
				print OUT "$$features{$id}{common}{groupID}\t$id\t$$features{$id}{common}{Description}\t$$features{$id}{common}{GN}\t$$features{$id}{common}{GPCRs}\t$$features{$id}{common}{TFs}\t$$features{$id}{common}{epigenetic_regularos}\t$$features{$id}{common}{kinases}\t$$features{$id}{common}{oncogenes}";
			} else {
				print OUT "$id\t$$features{$id}{common}{ModSites}\t$$features{$id}{common}{groupID}\t$$features{$id}{common}{proteinID}\t$$features{$id}{common}{Description}\t$$features{$id}{common}{GN}";
			}
		} elsif (defined($output_mode) and $output_mode==2 and $wNA==0) {
			print OUT "$id";
		}

		# batch specific features
	#	foreach my $batch (sort {$a <=> $b} keys %{$batches}) {
		my @tmp; foreach (keys %{$batches}) { s/batch//; push @tmp, $_; }
		foreach (sort {$a <=> $b} @tmp) {
			my $batch = "batch$_";
			if (defined($$features{$id}{batchSpecific}{$batch}{PSMs}) and 
			$$features{$id}{batchSpecific}{$batch}{PSMs}>0) { # quan. in this batch
				$bk++;
				push @IDbatches,$batch;

				# print annotation columns
				if (!defined($output_mode) or $output_mode==0 or $output_mode==1) {
					if ($input_mode==1) {
						print OUT "\t$batch\t$$features{$id}{batchSpecific}{$batch}{PSMs}\t$$features{$id}{batchSpecific}{$batch}{pepT}\t$$features{$id}{batchSpecific}{$batch}{pepU}\t$$features{$id}{batchSpecific}{$batch}{Coverage}";
					} else {
						print OUT "\t$batch\t$$features{$id}{batchSpecific}{$batch}{PSMs}";
					}
				}

				# print abundance columns
				unless (defined($output_mode) and $output_mode==2 and $wNA==1) {
					#my $n=$$parahash{'batches'}{$batch}{n};
					my $n=scalar(@{$$batches{$batch}});
					for (my $i=0;$i<$n;$i++) {
						if ($Om{$id}{$batch}{$batches{$batch}[$i]}) {
							print OUT "\t$$expression{$id}{$batch}{$$batches{$batch}[$i]}";
						} else {
							print OUT "\tNA";
						}
					}
				}
				#foreach my $sample (keys %{$$expression{$id}{$batch}}) {
					#print OUT "\t$$expression{$id}{$batch}{$sample}";
				#}
			} else { # missed in this batch
				if (!defined($output_mode) or $output_mode==0 or $output_mode==1) {
					if ($input_mode==1) {
						print OUT "\tNA\t0\t0\t0\t0";
					} else {
						print OUT "\tNA\t0";
					}
					my $N=scalar(@{$$batches{$batch}});
					print OUT "\tNA"x$N;
				}
			}
		}
		
		unless (defined($output_mode) and $output_mode==2 and $wNA==1) {
			print OUT "\n";
		}

		# record in how many batches this protein is ID
		if (defined($counts{$bk})) {
			$counts{$bk}++;
		} else {
			$counts{$bk}=1;
		}

		# for batch specific protein, count numbers
		if ($bk==1) {
			$singleBatchCounts{shift @IDbatches}++;
		}
	}

	close OUT;

	# print header 
	if (!defined($output_mode) or $output_mode==0 or $output_mode==1) {
		my $annotationSign='';
		if (defined($output_mode) and $output_mode==1) {
			$annotationSign='# ';
		}

		open (OUT,">.tmp_header");
		print OUT $annotationSign,"Total protein number combining all batches: $totalCount\n";
		foreach my $bk (sort {$b <=> $a} keys %counts) {
			print OUT $annotationSign,"  In any $bk batche(s): $counts{$bk}\n";
		}
		print OUT $annotationSign,"Batch-specific protein numbers:\n";
	#	foreach my $batch (sort {$a <=> $b} keys %singleBatchCounts) {
	my @tmp; foreach (keys %singleBatchCounts) { s/batch//; push @tmp, $_; }
	foreach (sort {$a <=> $b} @tmp) {
		my $batch = "batch$_";
			print OUT $annotationSign,"  $batch: $singleBatchCounts{$batch}\n";
		}
		print OUT $annotationSign,"\n";
		close OUT;

		# combine
		system(qq(cat .tmp_header .tmp_body > $output; rm .tmp_header .tmp_body));
	} elsif (defined($output_mode) and $output_mode==2) {
		system(qq(mv .tmp_body $output));
	}
}

sub featureInitialization {
	my ($ft2name)=@_;
	$$ft2name{common}{groupID}='Protein Group#';
	$$ft2name{common}{Description}='Protein Description';
	$$ft2name{common}{proteinID}='Protein Accession #';
	$$ft2name{common}{GN}='GN';
	$$ft2name{common}{GPCRs}='GPCRs';
	$$ft2name{common}{TFs}='TFs';
	$$ft2name{common}{ModSites}='Mod sites';
	#$$ft2name{common}{epigenetic_regularos}='epigenetic_regularos';
	$$ft2name{common}{epigenetic_regularos}='epigenetic_factors';
	$$ft2name{common}{kinases}='kinases';
	$$ft2name{common}{oncogenes}='oncogenes';
	$$ft2name{batchSpecific}{PSMs}='PSM#';
	$$ft2name{batchSpecific}{pepT}='Total Peptide#';
	$$ft2name{batchSpecific}{pepU}='Unique Peptide#';
	$$ft2name{batchSpecific}{Coverage}='%Coverage ';
	#$$ft2name{batchSpecific}{}='Peptide of the Highest Score';
	#$$ft2name{batchSpecific}{}='';
}


sub parseHeader {
# build %batches{$batch}[0]='Sig126'
	my ($batch,$line,$parahash,$batches,$ft2name,$featureIndex)=@_;
	chomp($line);
	chop($line) if ($line =~ m/\r$/); # windows ending
	my @t=split /\t/,$line;
	my $n=$$parahash{'batches'}{$batch}{n};
	for (my $i=0;$i<$n;$i++) {
		$$batches{$batch}[$i]=$t[$#t-$n+1+$i];
	}

	# features
	for (my $i=0;$i<=$#t;$i++) {
		my $found=0; # whether matched already
		# common features
		foreach my $f (keys %{$$ft2name{common}}) {
			if ( $$ft2name{common}{$f} eq $t[$i] ) {
				$$featureIndex{common}{$f}=$i;
				$found=1;
				last;
			}
		}
		# batchSpecific
		unless ($found) {
			foreach my $f (keys %{$$ft2name{batchSpecific}}) {
				if ( $$ft2name{batchSpecific}{$f} eq $t[$i] ) {
					$$featureIndex{batchSpecific}{$f}=$i;
					last;
				}
			}
		}
	}
}

sub parse_params
{
	#my $TMT=10;
        my ($par,$parahash)=@_;
	#my $general=1;
	#my $batch='';
	my $normalization_standard=1;
        open(IN,$par);
        while (<IN>)
        {
                s/^\s+//;
                next if (/^#/);
                chomp;
		chop if (/\r$/); # windows ending

		s/\s*([;\#].*)$//; # delete comments
		next unless (/ = /);
		my ($key,$value) = split(' = ',$_);

		if ($key =~ m/batch/) {
			if ($key =~ m/input_path_batch(\d+)/) {
				my $batch = "batch$1";
				$$parahash{'batches'}{$batch}{path}=$value;
				#$$parahash{'batches'}{$batch}{n}=$TMT;
			} elsif ($key =~ m/input_n_batch(\d+)/) {
				my $batch = "batch$1";
				$$parahash{'batches'}{$batch}{n}=$value;
			} elsif ($key =~ m/internal_standard_batch(\d+)/) {
				my $batch = "batch$1";
				if ($normalization_standard) {
					$$parahash{'normalization_standard'}=$batch;
					$normalization_standard=0;
				}
				$$parahash{'batches'}{$batch}{internal_standarad}=$value;
			}
		} else {
			$$parahash{$key}=$value;
		}

=head
		if (/^batch inofmration/) {
			$general=0;
			next;
		}

                s/\s*([;\#].*)$//; # delete comments
		if ($general) {
			next unless (/ = /);
			my ($key,$value) = split(' = ',$_);
			$$parahash{$key}=$value;
		} else {
			if (/^batch/) {
				s/\:$//;
				$batch=$_;
				next;
			} else {
				next unless (/ = /);
				my ($key,$value) = split(' = ',$_);
				$$parahash{'batches'}{$batch}{$key}=$value;
			}
		}
                if ($key =~ /^batch/)
                {
                        $$parahash{'batches'}{$key}=$value;
                }
                else
                {
                        $$parahash{$key}=$value;
                }
=cut
        }
        close IN;

	$$parahash{'focus_species'}='HUMAN';
}

sub log10 {
     my $n = shift;
     return log($n)/log(10);
}

sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}


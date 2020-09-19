#!/usr/bin/perl

package Utils::Parse;
use strict;
use warnings;
use File::Basename;
use Cwd;

sub new {
	my ($class) = @_;
	my $self = {};
	bless ($self,$class);
	return $self;
}

sub parseParams {
	shift @_;
	my ($paramFile, $params) = @_;
	open (IN, "<", $paramFile) or die "  Cannot open $paramFile\n";
	my $nComparisons = 0;
	while (<IN>) {
		next if (!/\=/ || /^\#/ || /^\s+\#/);
		$_ =~ s/\s+//g;
		$_ =~ s/\#.*//;
		my ($key, $value)  = split("\=", $_);
		if ($key =~ /(comparison_groups)_(.*)/) {
			$$params{'comparisons'}{$nComparisons} = $value;	## sample labels to be compared
			$$params{'comparisonNames'}{$nComparisons} = $2;	## name of comparison study
			$nComparisons++;
		} else {
			$$params{$key} = $value;
		}
	}

	## Check the validity of parameters
	if (!defined $$params{'idtxt'} || !-e $$params{'idtxt'}) {
		die "  Please put a correct path for 'idtxt'\n";
	}
	if (!defined $$params{'id_all_prot_quan'}) {
		print "  Please put a correct parameter for 'id_all_prot_quan'\n";
		print "  If you do not want to use this parameter, please set 'id_all_prot_quan = 0'\n";
		exit;
	}
	if (defined $$params{'id_all_prot_quan'} && $$params{'id_all_prot_quan'} ne 0 && !-e $$params{'id_all_prot_quan'}) {
		print "  Please put a correct path for 'id_all_prot_quan' file\n";
	}
	if (!defined $$params{'save_dir'}) {
		die "  Please put a correct directory name for 'save_dir'\n";
	}
	if (!defined $$params{'ppi_filter'} || $$params{'ppi_filter'} < 0 || $$params{'ppi_filter'} > 100) {
		print "  Please put a correct parameter for 'ppi_filter'\n";
		print "  'ppi_filter' parameter should be between 0 and 100 (percentage)\n";
		exit;
	}
	if (!defined $$params{'min_intensity_method'}) {
		die "  Please put a correct parameter for 'min_intensity_method'\n";
	} else {
		my $methods = $$params{'min_intensity_method'};
		$methods =~ s/\s+//g;
		my @methodsArray = split(/,/, $methods);
		foreach (@methodsArray) {
			## parameter should be either 0, 1, 2, 3 or 4
			if ($_ != 0 && $_ != 1 && $_ != 2 && $_ != 3 && $_ != 4) {
				die "  Please put a correct parameter for 'min_intensity_method'\n";
			}
		}
	}
	if (!defined $$params{'min_intensity_value'}) {
		die "  Please put a correct parameter for 'min_intensity_value'\n";
	} else {
		my @methods = split(/,/, $$params{'min_intensity_method'});
		my @values =  split(/,/, $$params{'min_intensity_value'});
		my $nMethods = scalar(@methods);
		my $nValues = scalar(@values);
		## Number of methods and values should be the same
		if ($nMethods != $nValues) {
			die "  Please put a correct parameter for 'min_intensity_value'\n";
		}
	}
	if (defined $$params{'min_intensity_method_1_2_psm'}) {
		my $methods = $$params{'min_intensity_method_1_2_psm'};
		$methods =~ s/\s+//g;
		my @methodsArray = split(/,/, $methods);
		foreach (@methodsArray) {
			## parameter should be either 1, 2, 3 or 4
			if ($_ != 0 && $_ != 1 && $_ != 2 && $_ != 3 && $_ != 4) {
				die "  Please put a correct parameter for 'min_intensity_method_1_2_psm'\n";
			}
		}		
	}
	if (defined $$params{'min_intensity_value_1_2_psm'}) {
		my @methods = split(/,/, $$params{'min_intensity_method_1_2_psm'});
		my @values = split(/,/, $$params{'min_intensity_value_1_2_psm'});
		my $nMethods = scalar(@methods);
		my $nValues = scalar(@values);
		## Number of methods and values should be the same
		if ($nMethods != $nValues) {
			die "  Please put a correct parameter for 'min_intensity_value_1_2_psm'\n";
		}		
	}
	if(!defined $$params{'tmt_reporters_used'}) {
		die "  Please put a correct parameter for 'tmt_reporters_used'\n";
	}
	if(!defined $$params{'tmt_peak_extraction_second_sd'}) {
		die "  Please put a correct parameter for 'tmt_peak_extraction_second_sd'\n";
	}
	if(!defined $$params{'tmt_peak_extraction_method'}) {
		die "  Please put a correct parameter for 'tmt_peak_extraction_method'\n";
	} elsif ($$params{'tmt_peak_extraction_method'} != 1 && $$params{'tmt_peak_extraction_method'} != 2) {
		die "  Please put a correct parameter for 'tmt_peak_extraction_method'\n";
	}
	if (!defined $$params{'impurity_correction'}) {
		die "  Please put a correct parameter for 'impurity_correction'\n";
	} else {
		if ($$params{'impurity_correction'} != 0 && $$params{'impurity_correction'} != 1) {
			die "  Please put a correct parameter for 'impurity_correction'\n";
		}
	}
	if ($$params{'impurity_correction'} == 1) {
		## when 'impurity_correction' = 1, impurity_matrix should be correctly defined
		if (!defined $$params{'impurity_matrix'} || !-e $$params{'impurity_matrix'}) {
			die "  Please put a correct parameter for 'impurity_matrix'\n";
		}
	}
	if (!defined $$params{'comparison_analysis'}) {
		die "  Please put a correct parameter for 'comparison_analysis'\n";
	} else {
		if ($$params{'comparison_analysis'} != 0 && $$params{'comparison_analysis'} != 1) {
			die "  Please put a correct parameter for 'comparison_analysis'\n";
		}
	}
	if ($$params{'comparison_analysis'} == 1 && !defined $$params{'comparisons'}) {
		die "  'comparison_groups_...' should be specified for the comparison analysis\n";
	}
	if (!defined $$params{'loading_bias_correction'}) {
		die "  Please put a correct parameter for 'loading_bias_correction'\n";
	} else {
		if ($$params{'loading_bias_correction'} != 0 && $$params{'loading_bias_correction'} != 1) {
			die "  Please put a correct parameter for 'loading_bias_correction'\n";
		}
	}
	if ($$params{'loading_bias_correction'} == 1) {
		if (!defined $$params{'loading_bias_correction_method'}) {
			die "  Please put a correct parameter for 'loading_bias_correction_method'\n";
		} else {
			if ($$params{'loading_bias_correction_method'} != 1 && $$params{'loading_bias_correction_method'} != 2) {
				die "  Please put a correct parameter for 'loading_bias_correction_method'\n";
			}
		}
	}
	if (!defined $$params{'SNratio_for_correction'}) {
		die "  Please put a correct parameter for 'SNratio_for_correction'\n";
	}
	if (!defined $$params{'percentage_trimmed'}) {
		die "  Please put a correct parameter for 'percentage_trimmed'\n";
	}
	if (!defined $$params{'interference_removal'}) {
		$$params{'interference_removal'} = 0;
	}
	
	## Hard-code some default parameters
	$$params{'quan_method'} = "MS2";
	$$params{'impurity_batch'} = "Batch2";
	$$params{'outlier_threshold'} = "Q95";
	$$params{'FDR_method'} = "BH";
	
	## Check whether parameters are well matched	
	my @reportersUsed = split(/;/, $$params{'tmt_reporters_used'});
	my $nReportersUsed = scalar(@reportersUsed);
	my %sample2reporter;
	foreach my $reporterUsed (@reportersUsed) {
		## Check whether the sample labels are assigned to the reporters used
		if (!defined $$params{$reporterUsed}) {
			die "  Please specify a sample label corresponding to $reporterUsed\n";
		} else {
			$sample2reporter{$$params{$reporterUsed}} = $reporterUsed;
		}
	}
	my $tmt = basename($$params{'impurity_matrix'});	## TMT8 or TMT10
	($tmt) = ($tmt =~ /TMT(\d+)\.ini/);	 
	if ($tmt != $nReportersUsed && $$params{'impurity_correction'} == 1) {
		print "  For the impurity correction, all TMT-reporters should be used\n";
		print "  If you want to process a part of reporters, please set 'impurity_correction' parameter to 0\n";
		exit;
	}	
	if ($$params{'comparison_analysis'} == 1) {		
		foreach (keys %{$$params{'comparisons'}}) {
			my @comparedSamples = split(/[:,]/, $$params{'comparisons'}{$_});			
			foreach my $comparedSample (@comparedSamples) {
				if (!defined $sample2reporter{$comparedSample}) {
					print "  $comparedSample is not defined in the sample labels\n";
					print "  Please check sample labels for TMT reporters\n";
					exit;
				}
			}
		}
	}
	close (IN);
}

sub parseIdtxt{
	shift @_;
	my ($idTxt, $outfilehash, $peptidehash, $lowPPIs, $params, $logfile) = @_;	
	my ($outfilenum, $total) = (0, 0);

	open (ID, "<", $idTxt);
	## Check the database line
	my ($line, $variable);
	if ($line = <ID>) {
		if ($line =~ m/^Database/) {
			push (@{$lowPPIs}, $line);
			$line  =~ s/\s+//g;
			($variable, $$params{'database'}) = split(/\=/, $line);
		} else {
			die "  Wrong ID.txt format: database line missed\n";
		}
	} else {
		die " Empty ID.txt!!!\n";
	}
	
	## Check the header line
	if ($line = <ID>) {
		push (@{$lowPPIs}, $line);
		chomp ($line);
		my @parts = split(/;/, $line);
		if ($parts[0] ne 'Peptide' or scalar(@parts) < 15) {
			die "  Wrong ID.txt format: header line check failed\n";
		}
	} else {
		die " Empty ID.txt!!!\n";
	}

	## Retrieve outfile lines
	my (%totalOut, %ppiFilterOut);
	while (<ID>) {
		chomp;
		$_ =~ s/\/\//\//g;	## Change all double slashes (//) to single slash (/)
		my @parts = split(/;/, $_);
		my ($origpeptide, $protein, $outfile, $mh, $calcmh, $ppm, $xcorr, $dcn, $group, $subgroup);
		my $ppi = 100;
		if (scalar(@parts) > 16) {
			($origpeptide, $protein, $outfile, $mh, $calcmh, $ppm, $xcorr, $dcn, $group, $subgroup,$ppi) = 
			($parts[0], $parts[1], $parts[2], $parts[3], $parts[4], $parts[5], $parts[6], $parts[7], $parts[10], $parts[11], $parts[16]);
		} else {
   		 	($origpeptide, $protein, $outfile, $mh, $calcmh, $xcorr, $dcn, $group, $subgroup,$ppi) = 
   		 	($parts[0], $parts[1], $parts[2], $parts[3], $parts[4], $parts[5], $parts[6], $parts[9], $parts[10], $parts[15]);
		}
		$totalOut{$outfile} = '';
		$ppi =~ s/\%//;
		if (defined($$params{'ppi_filter'}) && $ppi =~ /\d+/) {
			if ($ppi < $$params{'ppi_filter'}) {
				$ppiFilterOut{$outfile} = '';
				push (@{$lowPPIs}, $_);
				next;
			}
		}

		## Get information for the outfilehash and peptidehash		
		$outfile =~ s/([\w\d\_\-]+)(\.\d+\.)(\d+)(\.)(\d)(\..*out)\Z/$1$2$3$4$5$6/;
		my $shortout = $1.$2.$3.$4.$5.$6;
  		my ($file, $scannum, $charge) = ($1, $3, $5);
		$scannum =~ s/^0+//;
		my $fraction = $outfile; 
		$fraction =~ s/($file\.\d+)\/$shortout/$1/; 
		my $path = $fraction;
		$path =~ s/($file\.\d+)//; 
		my $mzxmlfilename = $fraction; 
		$mzxmlfilename =~ s/\.\d+$/\.mzXML/;
		$$params{'mzXML'}{$file} = $mzxmlfilename;
				
		# Save information to outfilehash and peptidehash
		print "\r  Gathering information from $outfile      ";
		$scannum =~ s/^0// if ($scannum =~ /^0/);
		my ($leftAA, $peptide, $rightAA) = split(/\./, $origpeptide);
		$total++;
		$outfilenum++ if (!defined($$outfilehash{$file}{$outfile}));
		$$outfilehash{$file}{$outfile}{'leftAA'} = $leftAA;
		$$outfilehash{$file}{$outfile}{'rightAA'} = $rightAA;
   		$$outfilehash{$file}{$outfile}{'peptide'} = $peptide;
   		$$outfilehash{$file}{$outfile}{'proteins'}{$protein}++;
   		$$outfilehash{$file}{$outfile}{'charge'} = $charge;
   		$$outfilehash{$file}{$outfile}{'mh+'} = $mh;
   		$$outfilehash{$file}{$outfile}{'pepmw'} = $calcmh;
   		$$outfilehash{$file}{$outfile}{'xcorr'} = $xcorr;
   		$$outfilehash{$file}{$outfile}{'dcn'} = $dcn;
   		$$outfilehash{$file}{$outfile}{'scan'} = $scannum;
   		$$peptidehash{$outfile}{$peptide}{'protein'} = $protein;	
  	}
	 print $logfile "  Finished gathering information from outfiles      \n";			
 	 close (ID);

 	 print "\n";
	 print "  There are ",scalar(keys %totalOut)," PSMs in total\n";
	 print "  Removed ",scalar(keys %ppiFilterOut)," PSMs due to low precursor peak intensity ($$params{'ppi_filter'}\%)\n";
 	 print "  Gathered information from $outfilenum PSMs\n";
	 print $logfile "  There are ",scalar(keys %totalOut)," PSMs in total\n";
	 print $logfile "  Removed ",scalar(keys %ppiFilterOut)," PSMs due to low precursor peak intensity ($$params{'ppi_filter'}\%)\n";
 	 print $logfile "  Gathered information from $outfilenum PSMs\n";
}

sub parseJumpfPublicationTable {
	shift @_;
	my ($inputFile) = @_;	## JUMP -f publication table
	my %idHash;
	my @idAccArray;
	my $idHeader;
	open (ID, "<", $inputFile) or die "Cannot open $inputFile\n";	
	if (basename($inputFile) =~ /pep/) {	## id_xxx_pep.txt
		while (<ID>) {
			chomp;
			if ($_ =~ /^Peptide/) {
				$idHeader = $_;
				next;
			}
			next if (!defined ((split(/\t/, $_))[1]));
			my @elems = split(/\t/, $_);
			my $pep = (split(/\./, $elems[0]))[1];
			$idHash{$pep}{$elems[2]} = \@elems;
			push (@idAccArray, [$pep, $pep]);			
		}
	} elsif (basename($inputFile) =~ /prot/) {	## id_xxx_prot.txt
		<ID>;	## skip one line
		$idHeader = <ID>;
		chomp ($idHeader);
		while (<ID>) {
			chomp;
			my @elems = split(/\t/, $_);
			$idHash{$elems[1]} = \@elems;
			push (@idAccArray, [@elems[0, 1]]);
		}		
	} elsif (basename($inputFile) =~ /site/) {	## id_xxx_site.txt
		$idHeader = <ID>;
		chomp ($idHeader);
		while (<ID>) {
			chomp;
			my @elems = split(/\t/, $_);
			$idHash{$elems[1]} = \@elems;
			push (@idAccArray, [@elems[0, 1]]);
		}
	}
	close (ID);
	return ($idHeader, \%idHash, \@idAccArray);
}

sub parseNormalizedReporterIntensity {
	shift @_;
	my ($params, $refFile) = @_;
	my $saveDir = getcwd()."/".$$params{'save_dir'};
	my $prefix = (split(/\//, $saveDir))[-1];
	my $suffix = "prot";
	if (basename($refFile) =~ /pep/) {
		$suffix = "pep";
	} elsif (basename($refFile) =~ /site/) {
		$suffix = "site";
	}
	my $inputFile = $saveDir."/"."norm_".$prefix."_".$suffix.".txt";
	my %normHash;
	open (NORM, "<", $inputFile) or die "Cannot open $inputFile\n";
	my $header = <NORM>;
	chomp ($header);
	my @normReporters = split(/;/, $header);
	@normReporters = @normReporters[1..$#normReporters];
	while (<NORM>) {
		chomp;
		my @elems = split(/;/, $_);
		my @intensity;
		for (my $i = 1; $i < scalar(@elems); $i++) {
			push (@intensity, 2 ** $elems[$i]);
		}
		$normHash{$elems[0]} = \@intensity;  
	}
	close (NORM);
	return (\@normReporters, \%normHash);
}

sub parseComparisonResult {
	shift @_;
	my ($params, $refFile) = @_;
	my $saveDir = getcwd()."/".$$params{'save_dir'};
	my $prefix = (split(/\//, $saveDir))[-1];
	my $suffix = "prot";
	if (basename($refFile) =~ /pep/) {
		$suffix = "pep";
	} elsif (basename($refFile) =~ /site/) {
		$suffix = "site";
	}
	my %compHash;	
	foreach my $key (keys %{$$params{'comparisonNames'}}) {
		my $comparisonName = $$params{'comparisonNames'}{$key};
		my $compFile = $saveDir."/".$comparisonName."_".$prefix."_".$suffix.".txt";
		open (COMP, "<", $compFile) or die "Cannot open $compFile\n";
		my $header = <COMP>;
		chomp ($header);
		my @headerArray = split(/\t/, $header);
		@headerArray = @headerArray[1..$#headerArray];		
		$compHash{'header'}{$comparisonName} = \@headerArray;
		while (<COMP>) {
			chomp;
			my @elems = split(/\t/, $_);
			my @stats = @elems[1..$#elems];
			$compHash{$elems[0]}{$comparisonName} = \@stats;
		}
		close (COMP);
	}
	return (\%compHash);
}

1;

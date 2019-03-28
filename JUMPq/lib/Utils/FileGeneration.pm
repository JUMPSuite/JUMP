#!/usr/bin/perl

package Utils::FileGeneration;
use strict;
use warnings;
use File::Basename;
use Cwd;
use Statistics::Distributions;

sub new {
	my ($class) = @_;
	my $self = {};
	bless ($self,$class);
	return $self;
}

sub generateRawTxtFiles {
	shift @_;
	my ($outfileHash, $peptideHash, $pep2site, $params, $logFile) = @_;
	my $idTxt = $$params{'idtxt'};
	my $pho = 0;
	if (basename($idTxt) eq "IDmod.txt") {
		$pho = 1;
	}
	my @reporters = split(/\;/, $$params{'tmt_reporters_used'});
	my @reporterMzs = getReporterMass(@reporters);
	my $nReporters = scalar(@reporters);

	##########################################
	## 1st stage of filtering PSMs by	##
	## 1. low PPI 				##
	## 2. zero intensity			##
	## 3. min_intensity_method/value filter	##
	##########################################
	## Intialization of filters
	my @filterMethods = split(/,/, $$params{'min_intensity_method'});
	my @filterCutoffs = split(/,/, $$params{'min_intensity_value'});
	my @filterNames;
	for (my $i = 0; $i < scalar(@filterMethods); $i++) {
		if ($filterMethods[$i] == 1) {
			push (@filterNames, "minimum");
		} elsif ($filterMethods[$i] == 2) {
			push (@filterNames, "maximum");
		} elsif ($filterMethods[$i] == 3) {
			push (@filterNames, "mean");
		} elsif ($filterMethods[$i] == 4) {
			push (@filterNames, "median");
		} elsif ($filterMethods[$i] == 0) {
			push (@filterNames, "no filter");
		}
	}	
	my (%filtered, %nonFiltered);
	my (%zeroFiltered, %intFiltered);	
	my $nTotPSMs = scalar(keys %$peptideHash);
	my $nPSMs = 0;
	foreach my $fraction (keys %$outfileHash) {
		foreach my $psm (keys %{$$outfileHash{$fraction}}) {
			$nPSMs++;
			print "\r    Processing $nPSMs PSMs out of $nTotPSMs for the minimum intensity-based filtering";
			## Check PSMs filtered by zero intensity at least one channel
			my $peptide = $$outfileHash{$fraction}{$psm}{'peptide'};
			my @intensity;
			foreach my $reporterMz (@reporterMzs) {
				## When any peak is not found because of no m/z corresponding to the reporter
				if (!defined $$peptideHash{$psm}{$peptide}{'reporter'}{$reporterMz}{'mass'} || $$peptideHash{$psm}{$peptide}{'reporter'}{$reporterMz}{'mass'} == 0) {
					$zeroFiltered{$psm}++;
					$filtered{$psm}++;
					last;
				}
				## When any peak is not defined because of undefined or zero intensity
				## Any scan containing zero peaks after impurity correction will be included to SCANZERO file by this condition
				if (!defined $$peptideHash{$psm}{$peptide}{'reporter'}{$reporterMz}{'intensity'} || $$peptideHash{$psm}{$peptide}{'reporter'}{$reporterMz}{'intensity'} == 0) {
					$zeroFiltered{$psm}++;
					$filtered{$psm}++;
					last;
				}
				## Check PSMs filtered by low intensity (by minimum_intensity_method/value parameters)
				push (@intensity, $$peptideHash{$psm}{$peptide}{'reporter'}{$reporterMz}{'intensity'});
			}
			if (scalar(@intensity) == $nReporters) {
				my ($isFiltered, $methodFiltered) = intensityFiltering(\@intensity, \@filterMethods, \@filterCutoffs, 0);
				if ($isFiltered == 1) {	## Filtered due to the low intensity
					$intFiltered{$methodFiltered}{$psm}++;
					$filtered{$psm}++;
				} else {
					$nonFiltered{$psm} = 1;
				}
			}
		}
	}
	print "\n";
	print $logFile "    Processed $nPSMs PSMs out of $nTotPSMs for the minimum intensity-based filtering\n";
	
	my %psm_1_2_filtered;
	my (@filterMethods_1_2_psm, @filterCutoffs_1_2_psm, @filterNames_1_2_psm);
	my $nTotProts;
	if ($pho == 0) {
		##################################################
		## 2nd stage of filtering PSMs by		##
		## min_intensity_method/value_1_2_psm filter	##
		## valid only for whole proteome analysis	##
		##################################################
		## Intialization of filters 
		@filterMethods_1_2_psm = split(/,/, $$params{'min_intensity_method_1_2_psm'});
		@filterCutoffs_1_2_psm = split(/,/, $$params{'min_intensity_value_1_2_psm'});		
		for (my $i = 0; $i < scalar(@filterMethods_1_2_psm); $i++) {
			if ($filterMethods_1_2_psm[$i] == 1) {
				push (@filterNames_1_2_psm, "minimum");
			} elsif ($filterMethods_1_2_psm[$i] == 2) {
				push (@filterNames_1_2_psm, "maximum");
			} elsif ($filterMethods_1_2_psm[$i] == 3) {
				push (@filterNames_1_2_psm, "mean");
			} elsif ($filterMethods_1_2_psm[$i] == 4) {
				push (@filterNames_1_2_psm, "median");
			} elsif ($filterMethods_1_2_psm[$i] == 0) {
				push (@filterNames_1_2_psm, "no filter");
			}
		}
		## Create a protein-to-outfile mapping hash
		my %prot2psm;
		foreach my $fraction (keys %$outfileHash) {
			foreach my $psm (keys %{$$outfileHash{$fraction}}) {
				if (defined $nonFiltered{$psm}) {
					foreach my $protein (keys %{$$outfileHash{$fraction}{$psm}{'proteins'}}) {
						push (@{$prot2psm{$protein}}, $psm);
					}
				}
			}
		}
		## Check PSMs mapped to a protein
		$nTotProts = scalar(keys %prot2psm);
		my $nProts = 0;		
		foreach my $protein (keys %prot2psm) {
			$nProts++;
			print ("\r    Processing $nProts proteins out of $nTotProts to filter out unreliable PSMs");
			my @psms = uniq(@{$prot2psm{$protein}});
			my $n = scalar(@psms);
			if ($n == 1) {
				## 1. Protein mapped by only one PSM
				##    If the PSM is filtered by the filter(s), it will not be used to quantify the protein
				my @intArray;
				my $peptide = (keys %{$$peptideHash{$psms[0]}})[0];
				foreach my $reporterMz (@reporterMzs) {
					push (@intArray, $$peptideHash{$psms[0]}{$peptide}{'reporter'}{$reporterMz}{'intensity'});
				}
				my ($isFiltered, $methodFiltered) = intensityFiltering(\@intArray, \@filterMethods_1_2_psm, \@filterCutoffs_1_2_psm, 0);
				if ($isFiltered == 1) {
					## If the PSM is filtered, it will be removed from %nonFilteredPSMs and added to %filteredPSMs
					delete $nonFiltered{$psms[0]};
					$psm_1_2_filtered{$methodFiltered}{$psms[0]}++;
					$filtered{$psms[0]}++;
				}
			} elsif ($n == 2) {
				## 2. Protein mapped by two PSMs
				##    2.1. Apply minimum intensity-based filter(s)
				##         If none of PSMs passes the filter(s), both PSMs will be removed from %nonFilteredPSMs
				##         If one of PSM passes the filter(s), the other will be removed from %nonFilteredPSMs 
				##         If all of PSMs pass the filter(s), go to the next step
				##         2.1.1. For each PSM, Check the variation (i.e. stdev) across all reporters (take log2)
				##                One with more variation will be removed from %nonFilteredPSMs
				##                If both PSMs have the same variation, then one with smaller mean intensity will be removed
				my (@intArray0, @intArray1);
				my $peptide0 = (keys %{$$peptideHash{$psms[0]}})[0];
				my $peptide1 = (keys %{$$peptideHash{$psms[1]}})[0];
				foreach my $reporterMz (@reporterMzs) {
					push (@intArray0, $$peptideHash{$psms[0]}{$peptide0}{'reporter'}{$reporterMz}{'intensity'});
					push (@intArray1, $$peptideHash{$psms[1]}{$peptide1}{'reporter'}{$reporterMz}{'intensity'});
				}
				my ($isFiltered0, $methodFiltered0) = intensityFiltering(\@intArray0, \@filterMethods_1_2_psm, \@filterCutoffs_1_2_psm, 0);
				my ($isFiltered1, $methodFiltered1) = intensityFiltering(\@intArray1, \@filterMethods_1_2_psm, \@filterCutoffs_1_2_psm, 0);
				if ($isFiltered0 == 1 && $isFiltered1 == 1) {
					## Both PSMs are filtered
					delete $nonFiltered{$psms[0]};
					delete $nonFiltered{$psms[1]};
					$psm_1_2_filtered{$methodFiltered0}{$psms[0]}++;
					$psm_1_2_filtered{$methodFiltered1}{$psms[1]}++;
					$filtered{$psms[0]}++;
					$filtered{$psms[1]}++;
				} elsif ($isFiltered0 == 0 && $isFiltered1 == 1) {
					delete $nonFiltered{$psms[1]};
					$psm_1_2_filtered{$methodFiltered1}{$psms[1]}++;
					$filtered{$psms[1]}++;
				} elsif ($isFiltered0 == 1 && $isFiltered1 == 0) {
					delete $nonFiltered{$psms[0]};
					$psm_1_2_filtered{$methodFiltered0}{$psms[0]}++;
					$filtered{$psms[0]}++;
				} else {
					@intArray0 = map {log($_) / log(2)} @intArray0;
					@intArray1 = map {log($_) / log(2)} @intArray1;
					my $sigma0 = stdev(\@intArray0);
					my $sigma1 = stdev(\@intArray1);
					if ($sigma1 > $sigma0) {
						delete $nonFiltered{$psms[1]};
						$psm_1_2_filtered{'postFilter'}{$psms[1]}++;
						$filtered{$psms[1]}++;
					} elsif ($sigma0 > $sigma1) {
						delete $nonFiltered{$psms[0]};
						$psm_1_2_filtered{'postFilter'}{$psms[0]}++;
						$filtered{$psms[0]}++;
					} else {
						my $mu0 = mean(\@intArray0);
						my $mu1 = mean(\@intArray1);
						if ($mu0 > $mu1) {
							delete $nonFiltered{$psms[1]};
							$psm_1_2_filtered{'postFilter'}{$psms[1]}++;
							$filtered{$psms[1]}++;
						} else {
							delete $nonFiltered{$psms[0]};
							$psm_1_2_filtered{'postFilter'}{$psms[0]}++;
							$filtered{$psms[0]}++;
						}
					}
				}
			}
		}
		print "\n";
		print $logFile "    Processed $nProts proteins out of $nTotProts to filter out unreliable PSMs\n";
	}

	##########################
	## Write to files	##
	##########################
	## Open two txt files and write headers
	my %prots;
	my %peps;
	my %sites;
	my $name = $$params{'save_dir'};
	my $tmt = (split(/\//, $$params{'impurity_matrix'}))[-1];
	my $scantext = $name; 
	my $scanzero = $name;
	$scantext .= "\/raw_$name\_psm_nonzero.txt";
	$scanzero .= "\/raw_$name\_psm_zero.txt";
	open(SCANTXT, ">", $scantext);
	open(SCANZERO, ">", $scanzero);	
	open(IDFILE, "<", $idTxt) || die "  Cannot open the file: $idTxt\n";
	my $database = <IDFILE>;
	print SCANTXT $database;
	print SCANZERO $database;
	my $header = <IDFILE>;
	chomp $header;
	my $nonzeroHeader = $header . ";RT;K_y1;R_y1;charge;endingAA;nLabels";
	my $zeroHeader = $header;	
	print SCANTXT $nonzeroHeader;
	print SCANZERO $zeroHeader;
	for (my $i = 0; $i < $nReporters; $i++) {
		my $mzName = "mz" . $reporters[$i];
		$mzName =~ s/sig//g;
		$mzName = $mzName . " (" . $$params{$reporters[$i]} . ")";
		print SCANTXT ";$mzName";
		print SCANZERO ";$mzName";
	}
	for (my $i = 0; $i < $nReporters; $i++) {
		my $intName = $reporters[$i] . " (" . $$params{$reporters[$i]} . ")";
		print SCANTXT ";$intName";
		print SCANZERO ";$intName";
	}
	print SCANTXT "\n";
	print SCANZERO "\n";
	
	## Write to SCANTXT and SCANZERO
	while(<IDFILE>) {
		my $line = $_;
		chomp $line;
		$line =~ s/\/\//\//g;		 
		my @data = split(/\;/, $line);
		my ($pep, $prot, $psm) = @data[0..2];
		$pep = (split(/\./, $pep))[1];
		if (!defined $$peptideHash{$psm}) {
			print SCANZERO "$line;";
			foreach my $reportMz (@reporterMzs) {
				print SCANZERO "0;0;";
			}
			print SCANZERO "\n";
		} elsif (defined $filtered{$psm}) {
			print SCANZERO "$line;";
			foreach my $reporterMz (@reporterMzs) {
				if (defined $$peptideHash{$psm}{$pep}{'reporter'}{$reporterMz}{'mass'}) {
					print SCANZERO $$peptideHash{$psm}{$pep}{'reporter'}{$reporterMz}{'mass'},"\;";
				} else {
					print SCANZERO "0;";
				}
			}
			foreach my $reporterMz (@reporterMzs) {			
				if (defined $$peptideHash{$psm}{$pep}{'reporter'}{$reporterMz}{'intensity'}) {			
					print SCANZERO $$peptideHash{$psm}{$pep}{'reporter'}{$reporterMz}{'intensity'},"\;";
				} else {                        
					print SCANZERO "0;";
				}
			}
			print SCANZERO "\n";
		} else {
			$line = $line . ";" . $$peptideHash{$psm}{$pep}{'rt'} . ";" .
					$$peptideHash{$psm}{$pep}{'Ky1'} . ";" . 
					$$peptideHash{$psm}{$pep}{'Ry1'} . ";" . 
					$$peptideHash{$psm}{$pep}{'charge'} . ";" .
					$$peptideHash{$psm}{$pep}{'endingAA'} . ";" .				 
					$$peptideHash{$psm}{$pep}{'nLabels'};
			$prots{$prot}++;
			$peps{$pep}++;
			foreach my $site (@{$$pep2site{$pep}}) {
				$sites{$site}++;
			}			
			print SCANTXT "$line;";
			foreach my $reporterMz (@reporterMzs) {
				print SCANTXT $$peptideHash{$psm}{$pep}{'reporter'}{$reporterMz}{'mass'},"\;";
			}
			foreach my $reporterMz (@reporterMzs) {			
				print SCANTXT $$peptideHash{$psm}{$pep}{'reporter'}{$reporterMz}{'intensity'},"\;";
			}
			print SCANTXT "\n";
		}
	}
	close(SCANTXT);
	close(SCANZERO);

	## Summarize the filtering results
	my $nPepsToBeProcessed = scalar(keys %peps);
	my $nProtsToBeProcessed = scalar(keys %prots);	
	my $nSitesToBeProcessed = scalar(keys %sites);	
	print "    Among $nTotPSMs PSMs,\n";
	print $logFile "    Among $nTotPSMs PSMs,\n";
	my $nTotFiltered = 0;
	my $nZeros = scalar(keys %zeroFiltered);
	print "    Removed $nZeros PSMs due to the zero intensity at least one channel\n";
	print $logFile "    Removed $nZeros PSMs due to the zero intensity at least one channel\n";
	$nTotFiltered = $nZeros;
	if (@filterMethods) {
		for (my $i = 0; $i < scalar(@filterMethods); $i++) {
			my $n = scalar(keys %{$intFiltered{$filterNames[$i]}});
			print "    Removed $n PSMs due to the minimum intensity less than $filterCutoffs[$i] ($filterNames[$i] of each PSM)\n";
			print $logFile "    Removed $n PSMs due to the minimum intensity less than $filterCutoffs[$i] ($filterNames[$i] of each PSM)\n";
			$nTotFiltered += $n;
		}
	}
	if (@filterMethods_1_2_psm) {
		my $nProtsFiltered = $nTotProts - $nProtsToBeProcessed;
		print "    Further filtering of 1 or 2 PSMs mapped to a protein\n";
		print $logFile "    Further filtering of 1 or 2 PSMs mapped to a protein\n";
		for (my $i = 0; $i < scalar(@filterMethods_1_2_psm); $i++) {
			my $n = scalar(keys %{$psm_1_2_filtered{$filterNames_1_2_psm[$i]}});
			print "    Removed $n PSMs due to the minimum intensity less than $filterCutoffs_1_2_psm[$i] ($filterNames_1_2_psm[$i]_1_2_psm of each PSM)\n";
			print $logFile "    Removed $n PSMs due to the minimum intensity less than $filterCutoffs_1_2_psm[$i] ($filterNames_1_2_psm[$i]_1_2_psm of each PSM)\n";
			$nTotFiltered += $n;
		}
		my $n = scalar(keys %{$psm_1_2_filtered{'postFilter'}});
		print "    Removed $n PSMs due to the larger variation than the other PSM mapped to the same protein\n";
		print "    Removed $nProtsFiltered proteins due to the filtered PSMs\n";
		print $logFile "    Removed $n PSMs due to the larger variation than the other PSM mapped to the same protein\n";
		print $logFile "    Removed $nProtsFiltered proteins due to the filtered PSMs\n";
		$nTotFiltered += $n;
	}	
	my $nPSMsToBeProcessed = $nTotPSMs - $nTotFiltered;
	
	print "    Hereafter, $nPSMsToBeProcessed unique PSMs\n";
	print "               $nPepsToBeProcessed unique peptides\n";
	if ($pho == 0) {
		print "               $nProtsToBeProcessed total proteins\n";
	} elsif ($pho == 1 && $$params{'jump_l'} == 1) {
		print "               $nSitesToBeProcessed total protein modification sites\n";
	}
	print "               will be used for subsequent analyses\n";
	print $logFile "    Hereafter, $nPSMsToBeProcessed unique PSMs\n";
	print $logFile "               $nPepsToBeProcessed unique peptides\n";
	if ($pho == 0) {
		print $logFile "               $nProtsToBeProcessed total proteins\n";
	} elsif ($pho == 1 && $$params{'jump_l'} == 1) {
		print $logFile "               $nSitesToBeProcessed total protein modification sites\n";
	}
	print $logFile "               will be used for subsequent analyses\n";
}

sub generateNormalizedFiles {
	shift @_;
	## Read 'raw_.._scan.txt' file and create hashes for peptide -> outfile and protein -> outfile
	my ($rawFile, $pep2site, $params, $logFile) = @_;
	my (%rawScan, %normPsmIntensity, %normPepIntensity, %normProtIntensity, %normSiteIntensity);
	my (%pep2psm, %prot2psm, %site2pep, %psm2pep);
	
	## Whole proteome or phosphoproteome analysis
	my $pho = 0;
	if (basename($$params{'idtxt'}) eq "IDmod.txt") {
		$pho = 1;
	}
	
	my $loadingFile = basename($rawFile);
	print "    Loading $loadingFile\n";
	print $logFile "    Loading $loadingFile\n";
	open (RAW, "<", $rawFile) || die "  Cannot open $rawFile\n";
	<RAW>;	# Skip the first lines
	my $rawHeader = <RAW>;
	chomp ($rawHeader);
	$rawScan{'header'} = $rawHeader;
	while (<RAW>) {
		chomp;
		my @elems = split(/;/, $_);
		my ($pep, $prot, $psm) = @elems[0..2];
		$pep = (split(/\./, $pep))[1];
		$rawScan{$pep}{$psm} = $_;
		push (@{$pep2psm{$pep}}, $psm);
		$psm2pep{$psm} = $pep;
		## Whole proteome analysis requires a %prot2psm hash for summarizing PSMs into proteins
		## Phosphoproteome analysis requries a %site2pep hash for summarizing peptides into sites
		if ($pho == 0) {
			push (@{$prot2psm{$prot}}, $psm);
		} else {
			if (defined $$pep2site{$pep}) {
				my @sites = uniq (@{$$pep2site{$pep}});
				foreach my $site (@sites) {
					push (@{$site2pep{$site}}, $pep);
				}
			}
		}
	}
	close (RAW);

	## Read 'norm_..._psm.txt' file and create a hash for normalized intensities of each outfile	
	my $normFile = basename($rawFile);
	$normFile =~ s/raw/norm/;
	$normFile =~ s/\_nonzero//;
	$normFile = dirname($rawFile)."/".$normFile;
	$loadingFile = basename($normFile);
	print "    Loading $loadingFile\n";
	print $logFile "    Loading $loadingFile\n";
	open (NORM, "<", $normFile) || die "  Cannot open $normFile\n";
	my $header = <NORM>;
	chomp ($header);
	my $nReporters = 0;
	while (<NORM>) {
		chomp;
		my @elems = split(/;/, $_);
		$nReporters = scalar(@elems) - 1;
		$normPsmIntensity{$elems[0]} = [@elems[1..$#elems]];		
	}
	close (NORM);

	## Summarize PSM-intensities to peptide-intensity
	## It is used for both whole proteome and phosphoproteome analysis
	## 1. For each PSM, log2-intensities are mean-centered across reporters
	## 1.1. e.g. for a peptide mapped by 10 PSMs with 10 reporters,
	##			 10 x 10 intensity matrix will be given, where each row (= PSM) is mean-centered
	## 2. Take the column-mean intensity for each peptide-level reporter intensity
	## 2.1. Outlier(s) is/are removed using Dixon q-test
	##      e.g. for sig126, 10 x 1 column vector is used to obtain peptide-level sig126 intensity
	## 2.2. At this step, 1 x 10 (mean-centered) intensity for the peptide can be obtained
	## 3. Obtain a 'representative (log2-) abundance of the peptide by taking a global mean of top 3 PSMs
	## 3.1. e.g. among 10 PSMs, choose 3 most abundant PSMs (= 3 x 10 log2-intensity matrix, not mean-centered)
	##           then, take a global mean of the 3 x 10 matrix for the representative abundance
	## 4. Multiply the representative abundance to the mean-centered peptide-level intensity
	
	## Create a file for peptide-level normalized intensity
	my $nPeps = 0;
	my $nTotPeps = scalar(keys %pep2psm);
	my $pepNormFile = basename($normFile);
	$pepNormFile =~ s/psm\.txt/pep\.txt/;
	$pepNormFile = dirname($normFile) . "/" . $pepNormFile;
	open (PEP, ">", $pepNormFile) || die "  Cannot open $pepNormFile\n";
	$header =~ s/PSM/Peptide/;
	print PEP "$header\n";
	my $psm2pepNormFile = basename($normFile);
	$psm2pepNormFile =~ s/psm\.txt/psm_to_pep\.txt/;
	$psm2pepNormFile = dirname($normFile) . "/" . $psm2pepNormFile;
	open (PSM2PEP, ">", $psm2pepNormFile);
	print PSM2PEP "$rawHeader\n";
	my $psm2pepOutlierNormFile = basename($normFile);
	$psm2pepOutlierNormFile =~ s/psm\.txt/psm_to_pep_outlier\.txt/;
	$psm2pepOutlierNormFile = dirname($normFile) . "/" . $psm2pepOutlierNormFile;
	open (PSM2PEP_OUTLIER, ">", $psm2pepOutlierNormFile);
	print PSM2PEP_OUTLIER "$rawHeader\n";
	
	foreach my $peptide (keys %pep2psm) {
		$nPeps++;
		print "\r    Processing $nPeps peptides out of $nTotPeps";		
		my @pepIntensity;
		my @pepOutliers;
		my @psmIntensity;
		my @psmMeanIntensity;		
		my @psms = uniq(@{$pep2psm{$peptide}});
		my %indHash;
		for (my $i = 0; $i < scalar(@psms); $i++) {
			my $psm = $psms[$i];
			for (my $j = 0; $j < $nReporters; $j++) {
				$psmIntensity[$i][$j] = $normPsmIntensity{$psm}[$j];
			}
			my $psmMean = mean(\@{$psmIntensity[$i]});
			push (@psmMeanIntensity, $psmMean);
			@{$psmIntensity[$i]} = map {$_ - $psmMean} @{$psmIntensity[$i]};	## Mean-centering of PSM intensities across reporters
			$indHash{$i} = 0;	## Initialization of the hash 
		}
	
		## Outlier removal and summarization by taking a mean value		
		for (my $j = 0; $j < $nReporters; $j++) {
			my @nonFiltered;
			for (my $i = 0; $i < scalar(@psms); $i++) {			
				push (@nonFiltered, $psmIntensity[$i][$j]);
			}
			my @filteredInd = outlierTest(\@nonFiltered, $$params{'outlier_threshold'});
			foreach (@filteredInd) {
				$indHash{$_}++;
			}
		}	
		my @indRetained = sort {$a <=> $b} (grep {$indHash{$_} == $nReporters} (keys %indHash));
		my @indRemoved = sort {$a <=> $b} (grep {$indHash{$_} != $nReporters} (keys %indHash));
		if (@indRemoved) {
			my @removedPsms = @psms[@indRemoved];
			foreach my $psm (@removedPsms) {		
				my @elems = split(/;/, $rawScan{$peptide}{$psm});
				splice(@elems, scalar(@elems) - $nReporters, $nReporters, @{$normPsmIntensity{$psm}});
				print PSM2PEP_OUTLIER join(";", @elems), "\n";
			}
		}		
		next if (!@indRetained);
		for (my $j = 0; $j < $nReporters; $j++) {
			my @filteredIntensity;
			foreach my $i (@indRetained) {
				push (@filteredIntensity, $psmIntensity[$i][$j]);
			}
			push (@pepIntensity, mean(\@filteredIntensity));
		}
		
		## Calculate representative abundance of the peptide
		@psmMeanIntensity = sort {$b <=> $a} @psmMeanIntensity;
		if (scalar(@psmMeanIntensity) > 3) {
			@psmMeanIntensity = @psmMeanIntensity[0..2];
		}
		my $pepRepAbundance = mean(\@psmMeanIntensity);
		
		## Restore peptide intensities to log2-scale
		@pepIntensity = map {$_ + $pepRepAbundance} @pepIntensity;
		$normPepIntensity{$peptide}{'normIntensity'} = \@pepIntensity;
		
		## Print to files
		print PEP $peptide, ";", join(";", @pepIntensity), "\n";
		@psms = @psms[@indRetained];
		foreach my $psm (@psms) {		
			my @elems = split(/;/, $rawScan{$peptide}{$psm});
			splice(@elems, scalar(@elems) - $nReporters, $nReporters, @{$normPsmIntensity{$psm}});
			print PSM2PEP join(";", @elems), "\n";
		}
	
	}
	close (PEP);
	close (PSM2PEP);
	close (PSM2PEP_OUTLIER);
	print "\n";
	print $logFile "    Processed $nPeps peptides out of $nTotPeps\n";

	## Create a file for protein- or protein modification site-level normalized intensity
	if ($pho == 0) {	## whole proteome analysis
		## Summarize PSM-intensities to a protein
		## The summarization scheme is the same as a peptide as above
		my $nProts = 0;
		my $nTotProts = scalar(keys %prot2psm);
		my $protNormFile = basename($normFile);
		$protNormFile =~ s/psm\.txt/prot\.txt/;
		$protNormFile = dirname($normFile)."/".$protNormFile;	
		open (PROT, ">", $protNormFile) || die "  Cannot open $protNormFile\n";
		$header =~ s/Peptide/Protein/;
		print PROT "$header\n";
		my $psm2protNormFile = basename($normFile);
		$psm2protNormFile =~ s/psm\.txt/psm_to_prot\.txt/;
		$psm2protNormFile = dirname($normFile) . "/" . $psm2protNormFile;
		open (PSM2PROT, ">", $psm2protNormFile);
		$rawHeader =~ s/Peptide;Protein/Protein;Peptide/;
		print PSM2PROT "$rawHeader\n";
		my $psm2protOutlierNormFile = basename($normFile);
		$psm2protOutlierNormFile =~ s/psm\.txt/psm_to_prot_outlier\.txt/;
		$psm2protOutlierNormFile = dirname($normFile) . "/" . $psm2protOutlierNormFile;
		open (PSM2PROT_OUTLIER, ">", $psm2protOutlierNormFile);		
		print PSM2PROT_OUTLIER "$rawHeader\n";
		
		foreach my $protein (keys %prot2psm) {
			$nProts++;
			print "\r    Processing $nProts proteins out of $nTotProts";
			my @protIntensity;
			my @protOutliers;
			my @psmIntensity;
			my @psmMeanIntensity;
			my @psms = uniq(@{$prot2psm{$protein}});
			my %indHash;
			for (my $i = 0; $i < scalar(@psms); $i++) {
				my $psm = $psms[$i];
				for (my $j = 0; $j < $nReporters; $j++) {
					$psmIntensity[$i][$j] = $normPsmIntensity{$psm}[$j];
				}
				my $psmMean = mean(\@{$psmIntensity[$i]});
				push (@psmMeanIntensity, $psmMean);
				@{$psmIntensity[$i]} = map {$_ - $psmMean} @{$psmIntensity[$i]};	## Mean-centering of PSM intensities across reporters
				$indHash{$i} = 0;	## Initialization of the hash 
			}
			
			## Outlier removal and summarization by taking a mean value			
			for (my $j = 0; $j < $nReporters; $j++) {
				my @nonFiltered;
				for (my $i = 0; $i < scalar(@psms); $i++) {
					push (@nonFiltered, $psmIntensity[$i][$j]);
				}
				my @filteredInd = outlierTest(\@nonFiltered, $$params{'outlier_threshold'});
				foreach (@filteredInd) {
					$indHash{$_}++;
				}
			}
			my @indRetained = sort {$a <=> $b} (grep {$indHash{$_} == $nReporters} (keys %indHash));
			my @indRemoved = sort {$a <=> $b} (grep {$indHash{$_} != $nReporters} (keys %indHash));
			if (@indRemoved) {
				my @removedPsms = @psms[@indRemoved];
				foreach my $psm (@removedPsms) {		
					my @elems = split(/;/, $rawScan{$psm2pep{$psm}}{$psm});
					splice(@elems, scalar(@elems) - $nReporters, $nReporters, @{$normPsmIntensity{$psm}});
					print PSM2PROT_OUTLIER "$protein;$elems[0];" . join(";", @elems[2..$#elems]), "\n";
				}
			}
						
			if (!@indRetained) {
				next;
			}
			for (my $j = 0; $j < $nReporters; $j++) {
				my @filteredIntensity;
				foreach my $i (@indRetained) {
					push (@filteredIntensity, $psmIntensity[$i][$j]);
				}
				push (@protIntensity, mean(\@filteredIntensity));
			}
			
			## Calculate representative abundance of the protein		
			@psmMeanIntensity = sort {$b <=> $a} @psmMeanIntensity;
			if (scalar(@psmMeanIntensity) > 3) {
				@psmMeanIntensity = @psmMeanIntensity[0..2];					
			}
			my $protRepAbundance = mean(\@psmMeanIntensity);
			
			# Restore protein intensities to log2-scale
			@protIntensity = map {$_ + $protRepAbundance} @protIntensity;
			$normProtIntensity{$protein}{'normIntensity'} = \@protIntensity;
			
			## Print to files
			print PROT $protein, ";", join(";", @protIntensity), "\n";
			@psms = @psms[@indRetained];
			foreach my $psm (@psms) {
				my @elems = split(/;/, $rawScan{$psm2pep{$psm}}{$psm});
				splice(@elems, scalar(@elems) - $nReporters, $nReporters, @{$normPsmIntensity{$psm}});
				print PSM2PROT "$protein;$elems[0];" . join(";", @elems[2..$#elems]), "\n";
			}
		}
		close (PROT);
		close (PSM2PROT);
		close (PSM2PROT_OUTLIER);
		print "\n";
		print $logFile "    Processed $nProts proteins out of $nTotProts\n";
	} elsif ($pho == 1 && $$params{'jump_l'} == 1) {	## phosphoproteome analysis
		## Summarize peptide-intensities to a site-intensity
		## 1. For each site, get peptides corresponding to the modification site (%pep2site hash)
		## 1.1. e.g. three peptides can define sp|Q12345|proteinX:S20
		##           KYLLS#PVLM
		##           KYLLS#PVLMS#VLLS
		##           KYLLS#PVLMSVLLS
		## 2. Add up the intensities of the peptides in raw-scale and then take log2 for the normalized intensity
		## 2.1. e.g.                          sig126    sig127N ... 
		##           KYLLS#PVLM               13.50     13.99 ...
		##           KYLLS#PVLMS#VLLS         15.92     16.04 ...
		##           KYLLS#PVLMSVLLS          14.38     15.22 ...
		##           -> convert to raw-scale for adding up
		##           KYLLS#PVLM               11585     16271 ...
		##           KYLLS#PVLMS#VLLS         62001     67378 ...
		##           KYLLS#PVLMSVLLS          21321     38166 ...
		##           -> site intensity
		##           sp|Q12345|proteinX:S20	  94907     121815 ...
		##           -> convert back to log2-scale
		##           sp|Q12345|proteinX:S20	  16.53     16.89 ...	
		my $nSites = 0;
		my $nQuantifiedSites = 0;
		my $nSkippedPeptides = 0;
		my $nTotSites = scalar(keys %site2pep);
		my $siteNormFile = basename($normFile);
		$siteNormFile =~ s/psm\.txt/site\.txt/;
		$siteNormFile = dirname($normFile)."/".$siteNormFile;	
		open (SITE, ">", $siteNormFile) || die "  Cannot open $siteNormFile\n";
		$header =~ s/Peptide/Site/;
		print SITE "$header\n";
		foreach my $site (keys %site2pep) {
			$nSites++;
			print "\r    Processing $nSites protein modification sites out of $nTotSites";
			my @siteIntensity;
			my @siteOutliers;
			my @peptides = uniq(@{$site2pep{$site}});
			for (my $i = 0; $i < scalar(@peptides); $i++) {
				if (!defined $normPepIntensity{$peptides[$i]}) {
					$nSkippedPeptides++;
					next;
				}
				for (my $j = 0; $j < $nReporters; $j++) {
					$siteIntensity[$j] += 2 ** $normPepIntensity{$peptides[$i]}{'normIntensity'}[$j];
				}
			}
			if (@siteIntensity) {
				$nQuantifiedSites++;
				for (my $j = 0; $j < $nReporters; $j++) {
					$siteIntensity[$j] = log($siteIntensity[$j]) / log(2);	## convert back to log-scale
					push (@siteOutliers, 0);
				}
				$normSiteIntensity{$site}{'normIntensity'} = \@siteIntensity;
				$normSiteIntensity{$site}{'numOutliers'} = \@siteOutliers;
				print SITE $site, ";", join(";", @siteIntensity), "\n";		
			}
		}
		close (SITE);
		print "\n";
		print "    Finally, $nQuantifiedSites sites are quantified (according to outlier removal, $nSkippedPeptides peptides are not considered)\n";
		print $logFile "    Processed $nSites protein modification sites out of $nTotSites\n";		
	}
	
	print "\n";
	print $logFile "\n";
	return (\%rawScan, \%normPsmIntensity, \%normPepIntensity, \%normProtIntensity, \%normSiteIntensity);
}

sub generateSiteTableFiles {
	shift @_;
	## Create id_uni/all_site.txt files
	my ($params, $siteInfo) = @_;
	my $uniProtFile = dirname($$params{'idtxt'}) . "/publications/id_uni_prot.txt";
	my $allProtFile = dirname($$params{'idtxt'}) . "/publications/id_all_prot.txt";	
	my $intermediateDir = getcwd() . "\/$$params{'save_dir'}" . "/intermediate/";
	my $uniSiteFile = $intermediateDir . "id_uni_site.txt";
	my $allSiteFile = $intermediateDir . "id_all_site.txt";
	
	## If the files do not exist, create the files
	## id_uni_site.txt
	open (INPUT, "<", $uniProtFile) or die "  Cannot open $uniProtFile";
	open (OUTPUT, ">", $uniSiteFile) or die "  Cannot open $uniSiteFile";
	my $header = <INPUT>;
	chomp($header);
	my @colNames = split(/\t/, <INPUT>);
	$colNames[1] = "Protein modification site";
	splice(@colNames, 2, 0, "Modification Score");
	print OUTPUT join("\t", @colNames);
	while (<INPUT>) {
		chomp;
		my @elems = split(/\t/, $_);
		my $prot = $elems[1];
		foreach my $site (sort {$a cmp $b} keys %{$$siteInfo{$prot}}) {
			my @values = @elems;
			$values[1] = $site;
			splice(@values, 2, 0, $$siteInfo{$prot}{$site});
			print OUTPUT join("\t", @values), "\n";
		}
	}
	close (OUTPUT);
	close (INPUT);
	
	## id_all_site.txt
	open (INPUT, "<", $allProtFile) or die "  Cannot open $allProtFile";
	open (OUTPUT, ">", $allSiteFile) or die "  Cannot open $allSiteFile";
	$header = <INPUT>;
	chomp($header);
	@colNames = split(/\t/, <INPUT>);
	$colNames[1] = "Protein modification site";
	splice(@colNames, 2, 0, "Modification Score");
	print OUTPUT join("\t", @colNames);
	while (<INPUT>) {
		chomp;
		my @elems = split(/\t/, $_);
		my $prot = $elems[1];
		foreach my $site (sort {$a cmp $b} keys %{$$siteInfo{$prot}}) {
			my @values = @elems;
			$values[1] = $site;
			splice(@values, 2, 0, $$siteInfo{$prot}{$site});
			print OUTPUT join("\t", @values), "\n";
		}
	}
	close (OUTPUT);
	close (INPUT);
}

sub getReporterMass {
	my @reporters = @_;
	my $nReporters = scalar(@reporters);
	my @reporterMasses;
	for (my $i = 0; $i < $nReporters; $i++) {
		if ($reporters[$i] eq "sig126") {
			$reporterMasses[$i] = 126.127726;			
		} elsif ($reporters[$i] eq "sig127" || $reporters[$i] eq "sig127N") {
			$reporterMasses[$i] = 127.124761;
		} elsif ($reporters[$i] eq "sig127C") {
			$reporterMasses[$i] = 127.131081;
		} elsif ($reporters[$i] eq "sig128N") {
			$reporterMasses[$i] = 128.128116;
		} elsif ($reporters[$i] eq "sig128" || $reporters[$i] eq "sig128C") {
			$reporterMasses[$i] = 128.134436;
		} elsif ($reporters[$i] eq "sig129" || $reporters[$i] eq "sig129N") {
			$reporterMasses[$i] = 129.131471;
		} elsif ($reporters[$i] eq "sig129C") {
			$reporterMasses[$i] = 129.137790;			
		} elsif ($reporters[$i] eq "sig130N") {
			$reporterMasses[$i] = 130.134825;
		} elsif ($reporters[$i] eq "sig130" || $reporters[$i] eq "sig130C") {
			$reporterMasses[$i] = 130.141145;
		} elsif ($reporters[$i] eq "sig131" || $reporters[$i] eq "sig131N") {
			$reporterMasses[$i] = 131.138180;
		} elsif ($reporters[$i] eq "sig131C") {
			$reporterMasses[$i] = 131.1445001;
		} else {
			die "  $reporters[$i] is an incorrect reporter name\n";
		}	
	}
	return (@reporterMasses);
}

sub outlierTest {
	my ($input, $alpha) = @_;
	$alpha =~ s/[A-Za-z]//;
	my $n = scalar(@$input);
	my @ind;
	
	## Assume the maximum number of outliers is around 30% of number of samples
	my $nOutliers = int($n * 0.3 + 0.5);
	if ($nOutliers > $n - 2) {
		$nOutliers = $n - 2;	## if $nOutliers > $n - 2, esdCriticalValues cannot be calculated
	}
	if ($nOutliers > 1) {
		@ind = genESD($input, $alpha, $nOutliers);
	} else {
		if ($n > 10) {
			@ind = genESD($input, $alpha, $nOutliers);
		} else {
			@ind = qtest($input, $alpha);
		}
	}
	return (@ind);
}

sub qtest {
	my ($input, $alpha) = @_;
	my $Q_threshold = "Q".$alpha;
	my %Q_table = ('Q90' => {3 => 0.941, 4 => 0.765, 5 => 0.642, 6 => 0.560, 7 => 0.507, 8 => 0.468, 9 => 0.437, 10 => 0.412},
				   'Q95' => {3 => 0.970, 4 => 0.829, 5 => 0.710, 6 => 0.625, 7 => 0.568, 8 => 0.526, 9 => 0.493, 10 => 0.466},
				   'Q99' => {3 => 0.994, 4 => 0.926, 5 => 0.821, 6 => 0.740, 7 => 0.680, 8 => 0.634, 9 => 0.598, 10 => 0.568});
	die "No confidence found in Q_table\n" unless (exists $Q_table{$Q_threshold});
		
	## Remove any undefined value
	my @temp;
	foreach (@$input) {
		if (defined($_) && $_ ne "N/A") {	
			push (@temp, $_);
		}
	}
	my @input = @temp;	
	undef @temp;
	
	## Sort the input and shift, if necessary
	my @ind = sort {$input[$a] <=> $input[$b]} (0..$#input);	
	@input = @input[@ind];
	my $minValue = min(\@input);
	my @data;
	if ($minValue < 0) {
		@data = map {$_ - $minValue} @input;
	} else {
		@data = @input;
	}
	my $n = scalar(@data);
	if ($n >= 3) {
		my $range = $data[-1] - $data[0];
		my $gap1 = $data[1] - $data[0];
		my $gapN = $data[-1] - $data[-2];
		my $Qstat1 = $gap1 / $range;
		my $QstatN = $gapN / $range;
		if (!defined $Q_table{$Q_threshold}{$n} || !defined $Qstat1 || !defined $QstatN) {
			print;
		}
		if ($Qstat1 > $Q_table{$Q_threshold}{$n}) {
			splice(@ind, 0, 1);
		}
		if ($QstatN > $Q_table{$Q_threshold}{$n}) {
			splice(@ind, -1)
		}
	}
	return (@ind);
}

sub genESD {
	my ($data, $alpha, $ub) = @_;	
	my @y = @$data;
	my @y2 = @y;
	my $n = scalar(@y2);
	my $nOutliers = 0;	
	
	## Compute test statistics until up to r = $ub values have been removed from the data
	my @indRemoved;
	for (my $i = 0; $i < $ub; $i++) {
		my %RtoY2;
		my @Rs;
		my $mu = mean(\@y2);
		my $sigma = stdev(\@y2);
		for (my $j = 0; $j < scalar(@y2); $j++) {
			my $Rval = abs($y2[$j] - $mu) / $sigma;
			push (@Rs, $Rval);
			$RtoY2{$Rval} = $y2[$j];
		}
		my $Ri = max(\@Rs);
		
		## Index manipulation
		## Which value in the original data vector (@y) is removed?
		for (my $j = 0; $j < scalar(@y); $j++) {
			if ($y[$j] == $RtoY2{$Ri}) {
				push (@indRemoved, $j);
			}
		}
		## Reduction of @y2
		delete $RtoY2{$Ri};
		@y2 = values %RtoY2;
		
		## Compute critical value,
		my $critVal = esdCritical($alpha, $n, $i);
		if ($Ri > $critVal) {
			$nOutliers = $i + 1;
		}
	}	
	@indRemoved = @indRemoved[0..($nOutliers - 1)];
	my @indRetained;	## index to be retained in y
	if (@indRemoved) {
		my %indRemovedHash;
		foreach (@indRemoved) {
			$indRemovedHash{$_} = 1;
		}
		for (my $i = 0; $i < scalar(@y); $i++) {
			if (!defined $indRemovedHash{$i}) {
				push (@indRetained, $i);
			}
		}
	} else {
		@indRetained = (0..$#y);
	}
	return (@indRetained);
}

sub esdCritical {
	my ($alpha, $n, $i) = @_;
	$alpha = $alpha / 100;	## input alpha is the unit of percentage
	$alpha = 1 - $alpha;
	my $p = 1 - $alpha / (2 * ($n - $i));
	my $t = Statistics::Distributions::tdistr($n - $i - 2, 1 - $p);
	my $crit = ($t * ($n -$i -1)) / sqrt(($n - $i - 2 + $t ** 2) * ($n - $i));
	return ($crit);
}

sub min {
	my $data = shift;
	my $min = 10000000000;
	foreach (@$data) {
		$min = $_ if $_ < $min;
	}
	return $min;
}

sub max {
	my $data = shift;
	my $max = 0;
	foreach (@$data ) {
		$max = $_ if $_ > $max;
	}
	return $max;
}

sub median {
	my $data = shift;
    my @vals = sort {$a <=> $b} @$data;
    my $len = @vals;
    if ($len % 2) {
		return $vals[int($len / 2)];	# odd
    } else {
		return ($vals[int($len / 2)-1] + $vals[int($len / 2)]) / 2;	# even
    }
}

sub mean {
	my ($data) = @_;
	if (!@$data) {
		die "Empty array\n";
	}
	my $total = 0;
	my $n = 0;
	foreach (@$data) {
		$total += $_;
		$n++;
	}
	my $average = $total / $n;
	return $average;
}

sub stdev{
	my ($data) = @_;
	if (scalar(@{$data}) == 1) {
		return 0;
	}
	my $average = mean($data);
	my $sqtotal = 0;
	my $n = 0;
	foreach (@$data) {
		$sqtotal += ($average - $_) ** 2;
		$n++;
	}
	my $std = ($sqtotal / ($n - 1)) ** 0.5;
	return $std;
}

sub uniq {
    my %seen = ();
    grep { not $seen{$_}++ } @_;
}

sub intensityFiltering {
	my ($input, $methods, $cutoffs, $log2scale) = @_;
	my $nFilters = scalar(@$methods);
	my ($isFiltered, $methodFiltered) = (0, 0);
	my $minIntensity;
	for (my $i = 0; $i < $nFilters; $i++) {
		if ($$methods[$i] == 1) {
			$minIntensity = min(\@$input);
			if ($log2scale == 1) {
				$minIntensity = 2 ** $minIntensity;
			}
			if ($minIntensity < $$cutoffs[$i]) {
				$isFiltered = 1;
				$methodFiltered = "minimum";
				last;
			}
		} elsif ($$methods[$i] == 2) {
			$minIntensity = max(\@$input);
			if ($log2scale == 1) {
				$minIntensity = 2 ** $minIntensity;
			}
			if ($minIntensity < $$cutoffs[$i]) {
				$isFiltered = 1;
				$methodFiltered = "maximum";
				last;
			}
		} elsif ($$methods[$i] == 3) {
			$minIntensity = mean(\@$input);
			if ($log2scale == 1) {
				$minIntensity = 2 ** $minIntensity;
			}
			if ($minIntensity < $$cutoffs[$i]) {
				$isFiltered = 1;
				$methodFiltered = "mean";
				last;
			}
		} elsif ($$methods[$i] == 4) {
			$minIntensity = median(\@$input);
			if ($log2scale == 1) {
				$minIntensity = 2 ** $minIntensity;
			}
			if ($minIntensity < $$cutoffs[$i]) {
				$isFiltered = 1;
				$methodFiltered = "median";
				last;
			}
		} elsif ($$methods[$i] == 0) {
			last;
		} else {
			print "  minimum intensity-baed filtering failed\n";
			print "  please check the parameters related with minimum intensity-based filtering\n";
			exit;
		}
	}
	return ($isFiltered, $methodFiltered);
}

1;

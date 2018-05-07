#!/usr/bin/perl

package Utils::ReporterIonExtraction;
use strict;
use warnings;
use Cwd;
use File::Basename;

sub new {
	my ($class) = @_;
	my $self = {};
	bless ($self,$class);
	return $self;
}

sub getReporterIntensity {
	shift @_;
	my ($outfileHash, $ms2Hash, $peptideHash, $params, $logFile, $reporterSummary) = @_;
	my @reporters = split(/;/, $$params{'tmt_reporters_used'});
	my @reporterMzs = getReporterMass(@reporters);	
	foreach my $fraction (sort {$a cmp $b} keys %{$outfileHash}) {
		my ($nOutfiles, $nTotalOutfiles) = (0, scalar(keys %{$$outfileHash{$fraction}}));
		foreach my $outfile (keys %{$$outfileHash{$fraction}}) {
			## Extract TMT-reporter intensities
			$nOutfiles++;
			if (!defined($reporterSummary)) {
				print "\r  $fraction --> $nOutfiles of $nTotalOutfiles PSMs";
			}			
			my $scanNumber = $outfile;
			$scanNumber =~ s/(.*[\w\d\_\-]+)\.(\d+)\.(\d+)(\.)(\d)(\..*out)\Z/$2/;
			my $peptide = $$outfileHash{$fraction}{$outfile}{'peptide'};
			foreach my $reporterMz (@reporterMzs) {
				my $tol = 0;
				my $center = 0;
				if (!defined($reporterSummary)) {
					$center = $reporterMz;
					$tol = 10;	# hard-coded 10 ppm tolerance around the theoretical reporter mass				
				} else {
					$center = $reporterMz + $$reporterSummary{$reporterMz}{'average'} * $reporterMz / 1000000;
					$tol = $$reporterSummary{$reporterMz}{'std'} * $$params{'tmt_peak_extraction_second_sd'}; 
				}				
				my $uL = $center + $tol / 1000000 * $reporterMz;
				my $lL = $center - $tol / 1000000 * $reporterMz;
				my (@mzArray, @intArray);
				if (!defined $$ms2Hash{$fraction}{$scanNumber}) {
					$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} = 0;
					$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'intensity'} = 0;
					next;
				}
				for (my $i = 0; $i < scalar(@{$$ms2Hash{$fraction}{$scanNumber}{'mz'}}); $i++) {
					my $mz = $$ms2Hash{$fraction}{$scanNumber}{'mz'}[$i];
					my $int = $$ms2Hash{$fraction}{$scanNumber}{'int'}[$i];
					if ($mz < $lL) {
						next;
					} elsif ($mz >= $lL && $mz <= $uL) {
						push (@mzArray, $mz);
						push (@intArray, $int);
					} elsif ($mz > $uL) {
						last;
					}
				}
				if (scalar(@mzArray) == 0) {	## no peak matched to the reporter
					$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} = 0;
					$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'intensity'} = 0;
				} elsif (scalar(@mzArray) == 1) { ## only one peak matched to the reporter
					$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} = $mzArray[0];
					$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'intensity'} = $intArray[0];
				} else {	## multiple peaks matched to the reporters
					if ($$params{'tmt_peak_extraction_method'} == 1) {	## strongest peak
						my $ind = 0;
						($$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'intensity'}, $ind) = strongestPeakSelection(\@intArray);
						$$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} = $mzArray[$ind];
					} elsif ($$params{'tmt_peak_extraction_method'} == 2) {	## closest peak
						my $ind = 0;
						($$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'intensity'}, $ind) = closestPeakSelection(\@mzArray, $reporterMz, 0);
					} else {
						print "  Please specify the correct peak extraction method\n";
						exit;
					}
				}
			}
			
			## Extract K-TMT and R-y1 ion intensities
			$$peptideHash{$outfile}{$peptide}{'rt'} = $$ms2Hash{$fraction}{$scanNumber}{'rt'};
			my $tol = 10;
			$$peptideHash{$outfile}{$peptide}{'Ky1'} = 0;
			$$peptideHash{$outfile}{$peptide}{'Ry1'} = 0;
			for (my $i = 0; $i < scalar(@{$$ms2Hash{$fraction}{$scanNumber}{'mz'}}); $i++) {
				my $mz = $$ms2Hash{$fraction}{$scanNumber}{'mz'}[$i];
				my $int = $$ms2Hash{$fraction}{$scanNumber}{'int'}[$i];
				my $Kmz = 376.2757362992;	## m/z of K-TMT-y1 ion
				my $Rmz = 175.1189521741;	## m/z of R-y1 ion
				my $KuL = $Kmz + $tol / 1000000 * $Kmz;
				my $KlL = $Kmz - $tol / 1000000 * $Kmz;
				my $RuL = $Rmz + $tol / 1000000 * $Rmz;
				my $RlL = $Rmz - $tol / 1000000 * $Rmz;
				if ($mz < $RlL) {
					next;
				} elsif ($mz >= $RlL && $mz <= $RuL) {
					if ($int > $$peptideHash{$outfile}{$peptide}{'Ry1'}) {
						$$peptideHash{$outfile}{$peptide}{'Ry1'} = $int;
					}
				} elsif ($mz >= $KlL && $mz <= $KuL) {
					if ($int > $$peptideHash{$outfile}{$peptide}{'Ky1'}) {
						$$peptideHash{$outfile}{$peptide}{'Ky1'} = $int;
					}
				} elsif ($mz > $KuL) {
					last;
				}
			}
			
			## Extract charge state, ending AA and number of TMT-labeling sites
			my $charge = (split(/\./, basename($outfile)))[-2];
			my $endingAA = substr($peptide, -1);
			my $nLabels = ($peptide =~ tr/K//) + 1;
			$$peptideHash{$outfile}{$peptide}{'charge'} = $charge;
			$$peptideHash{$outfile}{$peptide}{'endingAA'} = $endingAA;
			$$peptideHash{$outfile}{$peptide}{'nLabels'} = $nLabels;
		}
		if (!defined($reporterSummary)) {
			print $logFile "  $fraction --> $nOutfiles of $nTotalOutfiles PSMs\n";
		}
	}		
	print "\n";
	print $logFile "\n";
}

sub strongestPeakSelection {
	my ($array) = @_;
	my ($val, $ind) = (0, 0);
	for (my $i = 0; $i < scalar(@{$array}); $i++) {
		if ($$array[$i] >= $val) {
			$val = $$array[$i];
			$ind = $i;
		}
	}
	return ($val, $ind);
}

sub closestPeakSelection {
	my ($array, $reporterMz, $massShift) = @_;
	my ($closestMz, $diff, $ind) = (0, 10000, 0);
	for(my $i = 0; $i < scalar(@{$array}); $i ++) {
		if (abs($$array[$i] - $reporterMz - $massShift) < $diff) {
			$ind = $i;
			$closestMz = $$array[$i];
			$diff = abs($$array[$i] - $reporterMz - $massShift);
		}
	}
	return ($closestMz, $ind);
}

sub refineReporterIntensity {
	my $self = shift @_;
	my ($outfileHash, $ms2Hash, $peptideHash, $params, $logFile) = @_;
	my %reporterSummary;
	my @reporters = split(/;/, $$params{'tmt_reporters_used'});
	my @reporterMzs = getReporterMass(@reporters);
	my $nReporters = scalar(@reporters);
		
	## Summary of reporter ion intensities
	print "  Reporter ions extraction summary for the 1st round (mass tolerance = 10 ppm)\n";
	print $logFile "  Reporter ions extraction summary for the 1st round (mass tolerance = 10 ppm)\n";
	for (my $i = 0; $i < $nReporters; $i++) {
		my $reporterMz = $reporterMzs[$i];
		my $n = 0;
		my $nTot = 0;
		foreach my $outfile (keys %{$peptideHash}) {
			foreach my $peptide (keys %{$$peptideHash{$outfile}}) {
				if (defined($$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'}) && $$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} != 0) {
					$n++;
				}
				$nTot++;
			}
		}		
		my $percentage = $n / $nTot * 100;
		printf("    %s \t%.0f (%.2f%%) matched\n", $reporters[$i], $n, $percentage);
		printf $logFile ("    %s \t%.0f (%.2f%%) matched\n",$reporters[$i], $n, $percentage);
	}
		
	my %massShiftHash = ();
	my $shiftPpm = 0;
	print "  Standard deviation of mass shift:\n";
	print $logFile "  Standard deviation of mass shift:\n";
	for (my $i = 0; $i < $nReporters; $i++) {
		my $reporterMz = $reporterMzs[$i];
		my @tempArray;
		foreach my $outfile (keys %{$peptideHash}) {
			foreach my $peptide (keys %{$$peptideHash{$outfile}}) {
				if (defined($$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'}) && $$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} != 0) {
					if ($$params{'quan_method'} eq "MS2" || $$params{'quan_method'} eq "ms2") {
						$shiftPpm = (($$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} - $reporterMz) / $reporterMz) * 1000000;
						push (@tempArray, $shiftPpm);
					} else {
						print "  Currently, MS2-based quantification is only available\n";
						print "  Please choose MS2-based quantification method\n";
						exit;
					}
				}
			}
		}
		## Update reporter intensities and %reporterSummary
		my $avg = 0;
		my $sd = 0;
		if (scalar(@tempArray) > 1) {
			$avg = mean(\@tempArray);
			$sd = stdev(\@tempArray);
		}
        $reporterSummary{$reporterMz}{'average'} = $avg;
        $reporterSummary{$reporterMz}{'std'} = $sd;		
        print "    $reporters[$i] \t";
        printf(" m/z shift=%0.4f (ppm)\t\tsd=%0.4f\n", abs($avg), $sd);
        print $logFile "    $reporters[$i] \t";
        printf $logFile (" m/z shift=%0.4f (ppm)\t\tsd=%0.4f\n", abs($avg), $sd);
	}	
	print "  Re-extracting reporter ions based on calibrated mass and defined SD mass tolerance\n";
	print $logFile "  Re-extracting reporter ions based on calibrated mass and defined SD mass tolerance\n";

	## Check SD for the 2nd round of reporter extraction
	## This procedure tries to prevent the interference of reporter peaks due to the large SD
	## For example, assume the case as follows
	## sig127N	m/z shift = 5ppm, sd = 10
	## sig127C	m/z shift = -10ppm, sd = 20
	## and params{'tmt_peak_extraction_second_sd'} = 8
	## A new window for sig127N 
	## = 127.124761 + [5 * 127.124761 / 1e6] +/- [(10 * 8) * 127.124761 / 1e6]
	## = 127.115227 ~ 127.135567
	## A new window for sig127C 
	## = 127.131081 + [-10 * 127.131081 / 1e6] +/- [(20 * 8) * 127.131081 / 1e6]
	## = 127.109469 ~ 127.150151
	## New windows for sig127N and sig127C are overlapped. So, the reporter peaks can be interfered with each other
	my $SDlimit = 1000;	## recommended SD
	for (my $i = 0; $i < $nReporters - 1; $i++) {
		my $reporterMz1 = $reporterMzs[$i];
		my $reporterMz2 = $reporterMzs[$i + 1];
		my $center1 = $reporterMz1 + $reporterSummary{$reporterMz1}{'average'} * $reporterMz1 / 1000000;
		my $center2 = $reporterMz2 + $reporterSummary{$reporterMz2}{'average'} * $reporterMz2 / 1000000;
		my $coeff1 = $reporterSummary{$reporterMz1}{'std'} * $reporterMz1 / 1000000;
		my $coeff2 = $reporterSummary{$reporterMz2}{'std'} * $reporterMz2 / 1000000;
		my $currSD = ($center2 - $center1) / ($coeff1 + $coeff2);
		if ($currSD < $SDlimit) {
			$SDlimit = $currSD;
		}
	}
	printf "  According to the mass-shifts of reporter peaks, 'tmt_peak_extraction_second_sd' should be smaller than %0.2f\n", $SDlimit;
	print "  Current 'tmt_peak_extraction_second_sd' is $$params{'tmt_peak_extraction_second_sd'}\n";
	printf $logFile "  According to the mass-shifts of reporter peaks, 'tmt_peak_extraction_second_sd' should be smaller than %0.2f\n", $SDlimit;
	print $logFile "  Current 'tmt_peak_extraction_second_sd' is $$params{'tmt_peak_extraction_second_sd'}\n";
	if ($$params{'tmt_peak_extraction_second_sd'} > $SDlimit) {
		print "  == WARNING ==\n";
		print "  Reporter peak extraction may be inaccurate due to the large value of 'tmt_peak_extraction_second_sd'\n";
		printf "  Please set the parameter to the value smaller than %0.2f\n", $SDlimit;
		print $logFile "  == WARNING ==\n";
		print $logFile "  Reporter peak extraction may be inaccurate due to the large value of 'tmt_peak_extraction_second_sd'\n";
		printf $logFile "  Please set the parameter to the value smaller than %0.2f\n", $SDlimit;
	} 
	
	## Refine reporter intensity
	$self -> getReporterIntensity($outfileHash, $ms2Hash, $peptideHash, $params, $logFile, \%reporterSummary);
	
	## Summary of reporter intensities after 2nd round of selection
	print "  Reporter ions extraction summary for the 2nd round compared to the 1st extraction\n";
	print $logFile "  Reporter ions extraction summary for the 2nd round compared to the 1st extraction\n";	
	for (my $i = 0; $i < $nReporters; $i++) {
		my $reporterMz = $reporterMzs[$i];
		my $n = 0;
		my $nTotal = 0;
		foreach my $outfile (keys %{$peptideHash}) {
			foreach my $peptide (keys %{$$peptideHash{$outfile}}) {
				if(defined($$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'}) && $$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'mass'} != 0) {
					$n++;
				}
				$nTotal++;
			}
		}		
		my $percentage = $n /$nTotal * 100;
		printf("    %s \t%.0f (%.2f%%) matched\n", $reporters[$i], $n, $percentage);
		printf $logFile ("    %s \t%.0f (%.2f%%) matched\n", $reporters[$i], $n, $percentage);		
	}
}

sub correctImpurity {
	shift @_;
	my ($peptideHash, $params, $logFile) = @_;
	
	## Define an impurity matrix for the correction
	my @reporters = split(/;/, $$params{'tmt_reporters_used'});
	my @reporterMzs = getReporterMass(@reporters);
	my $nReporters = scalar(@reporters);
	my $impurityFile;
	if (defined $$params{'impurity_matrix'}) {
		$impurityFile = $$params{'impurity_matrix'};
		print "  $impurityFile will be used for the correction of TMT reporter impurities\n";
		print $logFile "  $impurityFile will be used for the correction of TMT reporter impurities\n";
	} else {
		die "  Please put a correct parameter for 'impurity_matrix'\n";
	}
	
	## In case of TMT10, Batch2 is chosen by default
	## For TMT8 and TMT11, there's only one batch
	my $impurityBatch;
	if (basename($impurityFile) eq "TMT10.ini") {
		$impurityBatch = "Batch2";
	} else {
		$impurityBatch = 0;
	}
	my ($impurityMatrix, $n) = getPurityMatrix($impurityFile, $impurityBatch);
	if ($nReporters != $n) {
		print "  Number of TMT-reporters in the parameter file and $impurityFile does not match\n";
		print "  If you want to process a part of reporters, please set 'impurity_correction' parameter to 0\n";
		exit;
	}
	
	my $R = Statistics::R -> new();
	$R -> set('x', \@$impurityMatrix);
	$R -> set('y', $nReporters);
	my $nTotPSMs = scalar(keys %{$peptideHash});
	my $nPSMs = 0;
	foreach my $outfile (keys %{$peptideHash}) {
		foreach my $peptide (keys %{$$peptideHash{$outfile}}) {
			my (@mz, @intensity);
			foreach my $reporterMz (sort {$a <=> $b} keys %{$$peptideHash{$outfile}{$peptide}{'reporter'}}) {
				push (@intensity, $$peptideHash{$outfile}{$peptide}{'reporter'}{$reporterMz}{'intensity'});
				push (@mz, $reporterMz);
			}
			$nPSMs++;
			print  "\r  Impurity correction: processing $nPSMs PSMs out of $nTotPSMs";			
			$R -> set('I', \@intensity);
			$R -> run(q`x <- matrix(x, ncol = y)`,
					  q`xinv <- round(solve(t(x)), 5)`,
					  q`realint <- as.numeric(round(c(I) %*% t(xinv), 5))`);
			my $res = $R -> get('realint');
			## If the intensity is over-corrected (i.e. smaller than 50% of uncorrected intensity),
			## the corrected intensity will be bound to 50% of uncorrected intensity
			for (my $i = 0; $i < scalar(@{$res}); $i++) {
				if ($$res[$i] < 0.5 * $intensity[$i]) {
					$$res[$i] = 0.5 * $intensity[$i];
				}
				$$peptideHash{$outfile}{$peptide}{'reporter'}{$mz[$i]}{'intensity'} = $$res[$i];
			}
		}
	}	
	print "\n";
	print $logFile "  Impurity correction: processed $nPSMs PSMs out of $nTotPSMs\n";
}

sub getPurityMatrix {
	my ($iniFile, $batch) = @_;
	my $flag = 0;
	my @dataMatrix;	
	my $nReporters = 0;
	
	if ($batch eq "0") {	## TMT8.ini or TMT11.ini
		open (INI, "<", $iniFile) || die "  Cannot open the ini file: $iniFile\n";
		while(<INI>) {
			if ($_ =~ /PurityMatrix/) {
				$flag = 1;				
				next;
			}
			if ($flag == 1 && $_ !~ /\[/) {
				my @data = split(/\s+/, $_);
				$nReporters = 0;
				for (my $i = 1; $i < scalar(@data); $i++) {
					push (@dataMatrix, $data[$i]);
					$nReporters++;
				}
			} else {
				$flag = 0;
			}
		}
		close (INI);	
	} else {	## TMT8.ini or TMT11.ini
		my $batchName;
		my %dataHash;
		open (INI, "<", $iniFile) || die "  Cannot open the ini file: $iniFile\n";
		while(<INI>) {
			if ($_ =~ /PurityMatrix_(.*)\]/) {
				$flag = 1;
				$batchName = $1;
				next;
			}
			if ($flag == 1 && $_ !~ /\[/) {
				my @data = split(/\s+/, $_);
				$nReporters = 0;
				for (my $i = 1; $i < scalar(@data); $i++) {
					push (@{$dataHash{$batchName}}, $data[$i]);
					$nReporters++;
				}
			} else {
				$flag = 0;
			}
		}
		close (INI);
		@dataMatrix = @{$dataHash{$batch}};	
	}
	return (\@dataMatrix, $nReporters);
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
	foreach ( @$data ) {
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
		$sqtotal += ($average-$_) ** 2;
		$n++;
	}
	my $std = ($sqtotal / ($n - 1)) ** 0.5;
	return $std;
}

sub uniq {
    my %seen = ();
    grep { not $seen{$_}++ } @_;
}
1;

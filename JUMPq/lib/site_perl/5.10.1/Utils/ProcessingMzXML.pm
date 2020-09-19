#!/usr/bin/perl

package Utils::ProcessingMzXML;

use strict;
use warnings;
use Cwd;
use File::Basename;
use Utils::XMLParser;

sub new {
	my ($class) = @_;
	my $self = {};
	bless ($self,$class);
	return $self;
}

sub getMS3 {
	shift @_;
	my ($xmlFile, $scanNumberArray, $logFile) = @_;
	my $fileName = basename($xmlFile);	
	my $xml = new Utils::XMLParser();
	my ($XML, $indexArray) = readMzXML($xmlFile);
	
	my %quanHash;
	## Mapping from MS2 scans to MS3 scans (in terms of scan number)
	my %ms2PrecursorMzToScan;
	my %ms2To3;
	print "  Retrieving scan information from $fileName file and correcting mass shifts (please wait)\n";
	print $logFile "  Retrieving scan information from $fileName file and correcting mass shifts (please wait)\n";
	for (my $i = 0; $i < scalar(@$indexArray); $i++) {
		my $scanNum = $i + 1;
		my $msLevel = $xml -> get_MSLevel(*XML, $$indexArray[$i]);
		if ($msLevel == 1) {
			undef %ms2PrecursorMzToScan;
		} else {
			if ($msLevel == 2) {	## MS2 scan
				my ($precursorMz, $dummy1, $dummy2, $dummy3) = $xml -> get_PrecursorMZINTACT(*XML, $$indexArray[$i]);
				$ms2PrecursorMzToScan{$precursorMz} = $scanNum;
			} elsif ($msLevel == 0) {	## MS3 scan
				my ($precursorMz, $dummy1, $dummy2, $dummy3) = $xml -> get_PrecursorMZINTACT(*XML, $$indexArray[$i]);
				if (defined $ms2PrecursorMzToScan{$precursorMz}) {
					$ms2To3{$ms2PrecursorMzToScan{$precursorMz}} = $scanNum;
				}
			}
		}
	}
	foreach my $scanNumber (@$scanNumberArray) {
		my $ms3ScanNumber = $ms2To3{$scanNumber};	## MS3 scan number corresponding to MS2 scan number
		next if (!defined $ms3ScanNumber);
		my $index = $$indexArray[$ms3ScanNumber - 1];
		my @peaksArray;
		$xml -> get_Peaks(*XML, \@peaksArray, $index);
		my (@mzArray, @intArray);		
		for (my $i = 0; $i < scalar(@peaksArray) / 2; $i += 2) {
			last if ($peaksArray[$i] > 135);
			push (@mzArray, $peaksArray[$i]);
			push (@intArray, $peaksArray[$i + 1]);
		}		
		$quanHash{$scanNumber}{'mz'} = \@mzArray;
		$quanHash{$scanNumber}{'int'} = \@intArray;
	}	
	## Mass correction
	massCorrection(\%quanHash, $logFile);	
	return (\%quanHash);
}

sub getMS2 {
	shift @_;
	my ($xmlFile, $scanNumberArray, $logFile) = @_;
	my $fileName = basename($xmlFile);
	my %ms2Hash;
	my $xml = new Utils::XMLParser();
	my ($XML, $indexArray) = readMzXML($xmlFile);
	print "  Retrieving scan information from $fileName file and correcting mass shifts (please wait)\n";
	print $logFile "  Retrieving scan information from $fileName file and correcting mass shifts (please wait)\n";
	foreach my $scanNumber (@$scanNumberArray) {
		my $index = $$indexArray[$scanNumber - 1];
		$ms2Hash{$scanNumber}{'rt'} = $xml -> get_RT(*XML, $index);
		my @peaksArray;
		$xml -> get_Peaks(*XML, \@peaksArray, $index);		
		my (@mzArray, @intArray);		
		for (my $i = 0; $i < scalar(@peaksArray) / 2; $i += 2) {
			last if ($peaksArray[$i] > 460);
			push (@mzArray, $peaksArray[$i]);
			push (@intArray, $peaksArray[$i + 1]);
		}		
		$ms2Hash{$scanNumber}{'mz'} = \@mzArray;
		$ms2Hash{$scanNumber}{'int'} = \@intArray;		
	}	
	## Mass correction
	massCorrection(\%ms2Hash, $logFile);	
	return (\%ms2Hash);
}

sub readMzXML {
	my ($mzXML) = @_;
	open (XML, "<", $mzXML) or die "  Cannot open $mzXML\n";
	my $xml = new Utils::XMLParser();
	my $indexOffset = $xml->get_IndexOffset(*XML); 
	my ($indexArray, $lastScan) = $xml->get_IndexArray(*XML, $indexOffset);
	return (*XML, $indexArray);
}

sub massCorrection {
	my ($ms2Hash, $logFile) = @_;
	my $refMass = 126.1277259380;	## theoretical mass of sig126 
	my @massShifts;
	
	## Calculation of mass shifts between observed and theoretical masses (sig126)
	foreach my $scan (keys %$ms2Hash) {
		my $obsMass = getObservedMass($$ms2Hash{$scan}{'mz'}, $$ms2Hash{$scan}{'int'}, $refMass);
		if ($obsMass != 0) {
			my $shift = ($obsMass - $refMass) / $refMass * 1e6;
			push (@massShifts, $shift);	## unit of ppm
		}
	}	
	## Filter the calculated mass shifts
	## 1. More than 50% of spectra should be used to calculate mass shifts
	## 2. Remove 20% highest/lowest mass shifts to get a more reliable correction factor 
	my $nScans = scalar(keys %$ms2Hash);
	my $nMassShifts = scalar(@massShifts);
	my $measuredPercentage = $nMassShifts / $nScans;

	## More than 50% of spectra should be used to calculate mass shifts
	## Otherwise, mass-shift correction will not be performed
	if ($measuredPercentage >= 0.5) {
		@massShifts = sort {$a <=> $b} @massShifts;
		## Filter 20% of highest/lowest values
		my $nFiltered = int(0.2 * $nMassShifts);
		splice (@massShifts, 0, $nFiltered);	## Remove the lowest mass shifts
		splice (@massShifts, $nMassShifts - $nFiltered, $nFiltered);	## Remove the highest mass shifts		
	
		## Calculate a global "correction factor" by taking the mean value of
		## mass shifts calculated from MS2 scans (i.e. mean of @massShifts)
		my $meanMassShift = mean(\@massShifts);
		my $stdMassShift = stdev(\@massShifts);
		printf "    Calculated mass-shift: mean = %.4f ppm and SD = %.4f ppm\n", $meanMassShift, $stdMassShift; 
		printf $logFile "    Calculated mass-shift: mean = %.4f ppm and SD = %.4f ppm\n", $meanMassShift, $stdMassShift;
		my $correctionFactor = mean(\@massShifts);
		
		## Mass-shift correction of all MS2 spectra using the correction factor
		foreach my $scan (keys %$ms2Hash) {
			for (my $i = 0; $i < scalar(@{$$ms2Hash{$scan}{'mz'}}); $i++) {
				$$ms2Hash{$scan}{'mz'}[$i] = $$ms2Hash{$scan}{'mz'}[$i] / (1 + $correctionFactor / 1e6);
			}
		}
	} else {
		print "    WARNING\n";
		printf("    Only %.2f%% of spectra (%d of %d effective MS2) is used to calculate mass shifts\n", 
				$measuredPercentage * 100, $nMassShifts, $nScans);
		print "    Mass correction will not be performed when less than 50% of spectra is used\n";
	}	

}

sub getObservedMass {
	my ($mzArray, $intArray, $refMass) = @_;
	my $nScans = scalar(@$mzArray);
	my $lL = $refMass - $refMass * 50 / 1e6;
	my $uL = $refMass + $refMass * 50 / 1e6;
	my ($obsMass, $obsIntensity) = (0, 0);
	for (my $i = 0; $i < $nScans; $i++) {
		last if ($$mzArray[$i] > $uL);
		if ($$mzArray[$i] >= $lL && $$mzArray[$i] <= $uL) {
			if ($$intArray[$i] > $obsIntensity) {
				$obsMass = $$mzArray[$i];
				$obsIntensity = $$intArray[$i];
			}
		}		
	}
	return ($obsMass);
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

1;

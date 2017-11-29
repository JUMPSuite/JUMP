#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd;
use Statistics::R;
use Excel::Writer::XLSX;

## Custom libraries
use FindBin qw($Bin);
use lib "$Bin";
use Utils::FileGeneration;
use Utils::Parse;
use Utils::ProcessingMzXML;
use Utils::PublicationTable;
use Utils::ReporterIonExtraction;
use Utils::XMLParser;

##########################
## Initialization	##
##########################
my ($paramFile) = @ARGV;
my $Rlib = $Bin . "/R/";
my $cwd = getcwd();
if (!(-e $paramFile)){
	print "$cwd does not contain a valid parameter file\n";
	exit;
}
my $file = Utils::FileGeneration -> new();
my $parse = Utils::Parse -> new();
my $procXML = Utils::ProcessingMzXML -> new();
my $publication = Utils::PublicationTable -> new();
my $extract = Utils::ReporterIonExtraction -> new();

##################################
## Parse a parameter file	##
##################################
my %params;
$parse -> parseParams($paramFile, \%params);

## Create a save directory
## 1. Save directory: a parameter file, a log file and raw files
## 2. Intermediate directory: a loading-bias log file, norm_.. files and comparision files
## 3. Publications direcotry: publication tables
## 4. Simplified publication directory: selected publication tables
$params{'save_dir'} = 'quan_' . $params{'save_dir'};
my $saveDir = getcwd() . "\/$params{'save_dir'}";
system ("mkdir $saveDir > /dev/null 2>&1");
my $intermediateDir = $saveDir."/intermediate";
system ("mkdir $intermediateDir > /dev/null 2>&1");
my $publicationDir = $saveDir."/publications";
system ("mkdir $publicationDir > /dev/null 2>&1");
my $simplePublicationDir = $saveDir."/simplified_publications";
system ("mkdir $simplePublicationDir > /dev/null 2>&1");

## Create log files and put it to the save directory
my $logFile;
open($logFile, ">", "$saveDir/jump_q.log") or die " Cannot open a log file\n";
my $loadingBiasLogFile;
open ($loadingBiasLogFile, ">", "$intermediateDir/jump_q_loading_bias.log") or die "  Cannot open a loading-bias log file";

## Copy the parameter file to the folder
system("cp $paramFile $saveDir");

my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
print "  jump -q started: ";
printf "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
print $logFile "  jump -q started: ";
printf $logFile "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;

##################
## Parse ID.txt	##
##################
my $idTxt = $params{'idtxt'};
if (!-e($idTxt)) {
	print "  Please input the correct path of ID.txt file\n";
	exit;
}
print "\n";
print $logFile "\n";
if (!defined($params{'tmt_peak_extraction_method'})) {
	print "  Please specify the peak extraction method!\n\n";
	exit;
}

## Parse ID.txt (or IDmod.txt) file
my (%outfileHash, %peptideHash);
my @lowPPIs;
$parse -> parseIdtxt($idTxt, \%outfileHash, \%peptideHash, \@lowPPIs, \%params, $logFile);

## In case of phosphoproteome analysis, some additional files need to be parsed and created
## 1. ID.lscore file is required
## 2. id_uni_site.txt needs to be created based on id_uni_prot.txt and ID.lscore
## 3. id_all_site.txt needs to be created based on id_all_prot.txt and ID.lscore
## All those files are the result of jump -l
my $pho = 0;
if (basename($idTxt) eq "IDmod.txt") {
	$pho = 1;
}
## If jump -l results do not exist, only peptide-level quantification will be performed
## Check the existence of jump -l results and put it to a parameter
my (%siteInfo, %pep2site);
$params{'jump_l'} = 1;
if ($pho) {
	## Check whether the user is using jump -f result or jump -l result
	my $inputDirName = (split(/\//, dirname($idTxt)))[-1];
	if ($inputDirName =~ /^sum/) {
		print "\n  ==================== WARNING =================================\n";
		print "  jump -f result will be used to run jump -q\n";
		print "  Hereafter, only peptide-level quantification will be performed\n";
		print "  ==============================================================\n\n";
		print $logFile "\n  ==================== WARNING =================================\n";
		print $logFile "  jump -f result will be used to run jump -q\n";
		print $logFile "  Hereafter, only peptide-level quantification will be performed\n";
		print $logFile "  ==============================================================\n\n";
		sleep(2);
		$params{'jump_l'} = 0;
	} elsif ($inputDirName =~ /^loc/) {
		## Parse jump -l result (ID.lscore file	by default)
		## Since the jump -l result may have a different name from ID.lscore, jump_l.params file should be parsed
		my $command = "ls " . dirname($idTxt) . "/*.params";
		my $jumplParams = `$command`;
		chomp($jumplParams);
		my $inputFile;
		if (-e $jumplParams) {
			my $jumplResult;
			open (JUMPL, "<", $jumplParams);
			while (<JUMPL>) {
				chomp;
				if ($_ =~ /^Output/) {
					$_ =~ s/\s+//g;
					$_ =~ s/\#.*//;
					$jumplResult = (split(/=/, $_))[1];
					last;
				}
			}
			close (JUMPL);	
			$inputFile = dirname($idTxt) . "/" . $jumplResult;
		}		
		if (!-e $inputFile) {
			print "\n  ======================================= WARNING ================================================\n";
			print "  jump -l (localization) results are not found\n";
			print "  Hereafter, only peptide-level quantification will be performed\n\n";
			print "  If you already ran jump -l and want to include a protein modification site-level quantification,\n";
			print "  1. please use IDmod.txt file obtained from jump -l\n";
			print "  2. please make sure that there is jump -l parameter file in the jump -l result directory\n";
			print "  3. also, please make sure that there's only one parameter file in the directory\n";
			print "  ================================================================================================\n\n";
			print $logFile "\n  ====================================== WARNING =================================================\n";
			print $logFile "  jump -l (localization) results are not found\n";
			print $logFile "  Hereafter, only peptide-level quantification will be performed\n\n";
			print $logFile "  If you already ran jump -l and want to include a protein modification site-level quantification,\n";
			print $logFile "  1. please use IDmod.txt file obtained from jump -l\n";
			print $logFile "  2. please make sure that there is jump_l.params file in the jump -l result directory\n";
			print $logFile "  3. also, please make sure that there's only one parameter file in the directory\n";
			print $logFile "  ================================================================================================\n\n";
			sleep(2);
			$params{'jump_l'} = 0;
		} else {
			open (INPUT, "<", $inputFile);
			<INPUT>;
			$siteInfo{'header'} = <INPUT>;	## skip the first two lines
			while (<INPUT>) {
				chomp;
				my @elems = split(/;/, $_);
				my ($prot, $psm, $pep) = ($elems[1], $elems[2], $elems[$#elems - 1]);
				$pep = (split(/\./, $pep))[1];
				my @sites = split(/,/, $elems[$#elems]);
				foreach my $site (@sites) {
					my $score;
					($site, $score) = split(/:/, $site);
					my $protSite = $prot . ":" . $site;
					if (defined $siteInfo{$prot}{$protSite}) {
						if ($score >= $siteInfo{$prot}{$protSite}) {
							$siteInfo{$prot}{$protSite} = $score;
						}
					} else {
						$siteInfo{$prot}{$protSite} = $score;
					}
					push (@{$pep2site{$pep}}, $protSite);
				}
			}
			close (INPUT);
			
			## Check id_uni/all_site.txt files
			## If they do not exist, create them
			my $uniSiteFile = dirname($params{'idtxt'}) . "/publications/id_uni_site.txt";
			my $allSiteFile = dirname($params{'idtxt'}) . "/publications/id_all_site.txt";
			$file -> generateSiteTableFiles(\%params, \%siteInfo);
		}		
	} else {
		print "  jump -q input should be either jump -f or jump -l result\n";
		print "  please check your 'idtxt' parameter\n";
		exit;
	}
}

##########################################
## Retrieve reporter intensities	##
##########################################
my %ms2Hash;
foreach my $key (sort {$a cmp $b} keys %outfileHash) {
	## Extract scan numbers corresponding to %outfilehash
	my @scanNumberArray;
	foreach my $outfile (keys %{$outfileHash{$key}}) {
		$outfile = basename($outfile);
		push (@scanNumberArray, (split(/\./, $outfile))[1]);
	}
	@scanNumberArray = sort {$a <=> $b} @scanNumberArray;
	@scanNumberArray = uniq(@scanNumberArray);
	my $mzXML = $params{'mzXML'}{$key};
	$ms2Hash{$key} = $procXML -> getMS2($mzXML, \@scanNumberArray, $logFile);
}
print "\n  Extracting reporter ion intensities\n";
print $logFile "\n  Extracting reporter ion intensities\n";
$extract -> getReporterIntensity(\%outfileHash, \%ms2Hash, \%peptideHash, \%params, $logFile);
$extract -> refineReporterIntensity(\%outfileHash, \%ms2Hash, \%peptideHash, \%params, $logFile);
print "\n";
print $logFile "\n";
if ($params{'impurity_correction'} == 1) {
	print "  Correcting the isotopic impurity of reporter ions\n";
	$extract -> correctImpurity(\%peptideHash, \%params, $logFile);
} else {
	print "  No correction of the isotopic impurity of reporter ions\n";
}
print "\n";
print $logFile "\n";

##################################################################
## Generate raw_..._scan.txt and raw_..._scan_zero.txt files	##
##################################################################
print "  Examining extracted PSMs (it may take a while)\n";
print $logFile "  Examining extracted PSMs (it may take a while)\n";
$file -> generateRawTxtFiles(\%outfileHash, \%peptideHash, \%pep2site, \%params, $logFile);

##############################################################
## Loading-bias representation and correction (optional)	##
##############################################################
print "\n";
print $logFile "\n";

## Show loading-biases over reporter ions before normalization
print "  Analyzing the loading-biases BEFORE correction (normalization)\n";
print $logFile "  Analayzing the loading-biases BEFORE correction (normalization)\n";
my $rawScanTxt = $saveDir."/"."raw_".$params{'save_dir'}."_psm_nonzero.txt";
my $noiseLevel = 1000;	## Hard-coded noise level
my $R = Statistics::R -> new();
$R -> run(q`rm(list = ls())`);
$R -> set('inputFile', $rawScanTxt);
$R -> set('noiseLevel', $noiseLevel);
$R -> set('SNratio', $params{'SNratio_for_correction'});
$R -> set('pctRemoval', $params{'percentage_trimmed'});
my $showLoadingBias = $Rlib."showLoadingBias.R";
$R -> run_from_file($showLoadingBias);
my $reporters = $R -> get('reporters');
my $nRows = $R -> get('nRows');
my $meanIntensity = $R -> get('meanIntensity');
my $sdIntensity = $R -> get('sdIntensity');
my $semIntensity = $R -> get('semIntensity');
my $canShowLoadingBias = $R -> get('canShowLoadingBias');
$R -> stop;
if ($canShowLoadingBias == 0) {
	print "    Due to the limited number of available PSMs, loading-biases are neither calculated nor corrected\n";
	print "    In order to obtain/correct loading-biases, try to reduce \"SNratio\" parameter\n";
	print $logFile "    Due to the limited number of available PSMs, loading-biases are neither calculated nor corrected\n";
	print $logFile "    In order to obtain/correct loading-biases, try to reduce \"SNratio\" parameter\n";
	print $loadingBiasLogFile "    Due to the limited number of available PSMs, loading-biases are neither calculated nor corrected\n";
	print $loadingBiasLogFile "    In order to obtain/correct loading-biases, try to reduce \"SNratio\" parameter\n";
} else {
	printf "    Noise-level (preset) = %d\n", $noiseLevel;
	print "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
	print "    percentage_trimmed = $params{'percentage_trimmed'} (from the parameter file)\n";
	printf "    Therefore, intensities less than %d are ignored in the calculation of loading-biases\n", $noiseLevel * $params{'SNratio_for_correction'};
	print "    And, for each reporter, $params{'percentage_trimmed'}\% of most variable intensities are also ignored (i.e. trimmed)\n\n";
	print "  Loading biases based on trimmed-mean intensities are as follows\n";
	print "    Reporters\tMean[%]\tSD[%]\tSEM[%]\t\#PSMs\n";
	printf $logFile "    Noise-level (preset) = %d\n", $noiseLevel;
	print $logFile "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
	print $logFile "    percentage_trimmed = $params{'percentage_trimmed'} (from the parameter file)\n";
	printf $logFile "    Therefore, intensities less than %d are ignored in the calculation of loading-biases\n", $noiseLevel * $params{'SNratio_for_correction'};
	print $logFile "    And, for each reporter, $params{'percentage_trimmed'}\% of most variable intensities are also ignored (i.e. trimmed)\n\n";
	print $logFile "  Loading biases based on trimmed-mean intensities are as follows\n";
	print $logFile "    Reporters\tMean[%]\tSD[%]\tSEM[%]\t\#PSMs\n";
	printf $loadingBiasLogFile "    Noise-level (preset) = %d\n", $noiseLevel;
	print $loadingBiasLogFile "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
	printf $loadingBiasLogFile "    Therefore, intensities less than %d are ignored in the calculation of loading-biases\n", $noiseLevel * $params{'SNratio_for_correction'};
	print $loadingBiasLogFile "    And, for each reporter, $params{'percentage_trimmed'}\% of most variable intensities are also ignored (i.e. trimmed)\n\n";
	print $loadingBiasLogFile "  Loading biases based on trimmed-mean intensities are as follows\n";
	print $loadingBiasLogFile "    Reporters\tMean[%]\tSD[%]\tSEM[%]\t\#PSMs\n";
	for (my $i = 0; $i < scalar(@{$reporters}); $i++) {
		printf "    %s\t%.2f\t%.2f\t%.2f\t%d\n", $$reporters[$i], $$meanIntensity[$i], $$sdIntensity[$i], $$semIntensity[$i], $nRows;
		printf $logFile "    %s\t%.2f\t%.2f\t%.2f\t%d\n", $$reporters[$i], $$meanIntensity[$i], $$sdIntensity[$i], $$semIntensity[$i], $nRows;
		printf $loadingBiasLogFile "    %s\t%.2f\t%.2f\t%.2f\t%d\n", $$reporters[$i], $$meanIntensity[$i], $$sdIntensity[$i], $$semIntensity[$i], $nRows;
	}	
}
print "\n";
print $logFile "\n";
print $loadingBiasLogFile "\n";
close ($loadingBiasLogFile);

## Normalization (a.k.a. loading-bias correction, PSM-level)
my $doNormalization = 0;
my $methodNormalization = 0;
if ($params{'loading_bias_correction'} == 1 && $canShowLoadingBias == 1) {
	$doNormalization = 1;
	print "  Data is being normalized (i.e. loading-bias correction)\n";
	print $logFile "  Data is being normalized (i.e. loading-bias correction)\n";
	if ($params{'loading_bias_correction_method'} == 1) {	# trimmed-mean
		$methodNormalization = 1;
		print "  by trimmed-mean method which makes trimmed-mean intensities of reporters identical\n";
		printf "    Noise-level (preset) = %d\n", $noiseLevel;
		print "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
		print "    percentage_trimmed = $params{'percentage_trimmed'} (from the parameter file)\n";
		print "    Therefore, for the calculation of each reporter's trimmed-mean intensity,\n";
		printf "    intensities less than %d and", $noiseLevel * $params{'SNratio_for_correction'};
		print " $params{'percentage_trimmed'}\% of most variable intensities are ignored\n";
		print $logFile "  by trimmed-mean method which makes trimmed-mean intensities of reporters identical\n";
		printf $logFile "    Noise-level (i.e. minimum intensity) = %d\n", $noiseLevel;
		print $logFile "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
		print $logFile "    percentage_trimmed = $params{'percentage_trimmed'} (from the parameter file)\n";
		print $logFile "    Therefore, for the calculation of each reporter's trimmed-mean intensity,\n";
		printf $logFile "    intensities less than %d and", $noiseLevel * $params{'SNratio_for_correction'};
		print $logFile " $params{'percentage_trimmed'}\% of most variable intensities are ignored\n";
	} elsif ($params{'loading_bias_correction_method'} == 2) {	# trimmed-median
		$methodNormalization = 2;
		print "  by trimmed-median method which makes trimmed-median intensities of reporters identical\n";
		printf "    Noise-level (preset) = %d\n", $noiseLevel;
		print "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
		print "    percentage_trimmed = $params{'percentage_trimmed'} (from the parameter file)\n";
		print "    Therefore, for the calculation of each reporter's trimmed-median intensity,\n";
		printf "    intensities less than %d and", $noiseLevel * $params{'SNratio_for_correction'};
		print " $params{'percentage_trimmed'}\% of most variable intensities are ignored\n";
		print $logFile "  by trimmed-median method which makes trimmed-median intensities of reporters identical\n";
		printf $logFile "    Noise-level (i.e. minimum intensity) = %d\n", $noiseLevel;
		print $logFile "    signal-to-noise ratio (SNratio) = $params{'SNratio_for_correction'} (from the parameter file)\n";
		print $logFile "    percentage_trimmed = $params{'percentage_trimmed'} (from the parameter file)\n";
		print $logFile "    Therefore, for the calculation of each reporter's trimmed-median intensity,\n";
		printf $logFile "    intensities less than %d and", $noiseLevel * $params{'SNratio_for_correction'};
		print $logFile " $params{'percentage_trimmed'}\% of most variable intensities are ignored\n";
	} else {
		print "  'loading_bias_correction_method' should be properly specified\n";
		exit;
	}
} else {
	print "  Data is not normalized (i.e. no loading-bias correction)\n";
	print $logFile "  Data is not normalized (i.e. no loading-bias correction)\n";
}
if ($params{'interference_removal'} == 1) {
	print "\n";
	print "  Interference in TMT quantification is being removed (based on y1-ion intensities)\n";
}
$R = Statistics::R -> new();
$R -> run(q`rm(list = ls())`);
$R -> set('saveDir', $saveDir);
$R -> set('inputFile', $rawScanTxt);
$R -> set('doNormalization', $doNormalization);
$R -> set('methodNormalization', $methodNormalization);
$R -> set('noiseLevel', $noiseLevel);
$R -> set('SNratio', $params{'SNratio_for_correction'});
$R -> set('pctRemoval', $params{'percentage_trimmed'});
$R -> set('interference_removal', $params{'interference_removal'});
my $normalizationR = $Rlib."normalization.R";
$R -> run_from_file($normalizationR);
my $canRemoveInterference = $R -> get('canRemoveInterference');
my $corrVec = $R -> get('rho');
$R -> stop;
if ($params{'interference_removal'} == 1 && $canRemoveInterference == 0) {
	print "  Interference in TMT quantification is not removed due to the lack of informative PSMs\n";
}

## Correlation matrix of PSM-level data
print "\n";
print "  Correlation between reporters (based on log2 PSM-level intensities)\n";
print $logFile "\n";
print $logFile "  Correlation between reporters (based on log2 PSM-level intensities)\n";
printCorrMatrix($corrVec, \%params, $logFile);

##############################################################################
## Summarization of PSMs into peptides/proteins								##
## Generate files for the summarized intensities of peptides and proteins	##
##############################################################################
print "  Summarizing PSMs into peptides and proteins (or protein modification sites) (it may take a while)\n";
print $logFile "  Summarizing PSMs into peptides and proteins (or protein modification sites) (it may take a while)\n";
my ($rawScanInfo, $normPsmIntensity, $normPepIntensity, $normProtIntensity, $normSiteIntensity) = $file -> generateNormalizedFiles($rawScanTxt, \%pep2site, \%params, $logFile);

##########################
## Statistical testing	##
##########################
my @comparisons;
my @comparisonNames;
if (defined $params{'comparison_analysis'} && $params{'comparison_analysis'} == 0) {
	print "  No comparison analysis\n";
	print $logFile "  No comparison analysis\n";

} else {
	foreach my $compInd (keys %{$params{'comparisons'}}) {
		push (@comparisons, $params{'comparisons'}{$compInd});
		push (@comparisonNames, $params{'comparisonNames'}{$compInd});
	}
	
	## Running R script
	print "  Statistical testing(s) is(are) being performed (it may take a while)\n\n";
	print $logFile "  Statistical testing(s) is(are) being performed\n\n";
	$R = Statistics::R -> new();
	$R -> run(q`rm(list = ls())`);
	$R -> set('saveDir', $saveDir);
	$R -> set('suffix', $params{'save_dir'});
	$R -> set('fdrMethod', $params{'FDR_method'});
	$R -> set('comparisons', \@comparisons);
	$R -> set('comparisonNames', \@comparisonNames);
	my $statTestR = $Rlib."statTest.R";
	$R -> run_from_file($statTestR);
	$R -> stop;	
}

##########################
## Publication tables	##
##########################
## 1. Parse jump -f publication table for peptide / protein entries
## 2. Parse "norm_..._pep (or prot).txt" file for the normalized reporter intensities
## 3. Parse comparison results (e.g. twoGroups_..._pep (or prot).txt) for the statistical testing results
## 4. Generate publication tables (including NA-table and txt file)
my $rawScanZeroTxt = $rawScanTxt;
$rawScanZeroTxt =~ s/psm\_nonzero\.txt/psm\_zero\.txt/;

if ($pho == 1) {
	print "  Generating publication tables for unique peptides\n";
	print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	print $logFile "  Generating publication tables for unique peptides\n";
	print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	$publication -> generatePublicationTables(\%params, "id_uni_pep.txt", $rawScanZeroTxt, $publicationDir);
	print "  Generating publication tables for all peptides\n";
	print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	print $logFile "  Generating publication tables for all peptides\n";
	print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	$publication -> generatePublicationTables(\%params, "id_all_pep.txt", $rawScanZeroTxt, $publicationDir);
	if ($params{'jump_l'} == 1) {
		print "  Generating publication tables for the modification sites of unique proteins\n";
		print "    It may take several minutes when there are lots of PSMs, peptides and protein modification sites\n";
		print $logFile "  Generating publication tables for the modification sites of unique proteins\n";
		print $logFile "    It may take several minutes when there are lots of PSMs, peptides and protein modification sites\n";
		$publication -> generatePublicationTables(\%params, "id_uni_site.txt", $rawScanZeroTxt, $publicationDir);
		print "  Generating publication tables for the modification sites of all proteins\n";
		print "    It may take several minutes when there are lots of PSMs, peptides and protein modification sites\n";
		print $logFile "  Generating publication tables for the modification sites of all proteins\n";
		print $logFile "    It may take several minutes when there are lots of PSMs, peptides and protein modification sites\n";
		$publication -> generatePublicationTables(\%params, "id_all_site.txt", $rawScanZeroTxt, $publicationDir);		
	}
	if (defined $params{'id_all_prot_quan'} && $params{'id_all_prot_quan'} ne 0) {
		my $protFile = $params{'id_all_prot_quan'};
		print "  Generating a combined publication table for unique peptides and proteins\n";
		print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
		print $logFile "  Generating a combined publication table for unique peptides and proteins\n";
		print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
		$publication -> generateCombinedPublicationTable(\%params, $publicationDir, $protFile);
	}
} else {
	print "  Generating publication tables for unique peptides\n";
	print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	print $logFile "  Generating publication tables for unique peptides\n";
	print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	$publication -> generatePublicationTables(\%params, "id_uni_pep.txt", $rawScanZeroTxt, $publicationDir);
	print "  Generating publication tables for all peptides\n";
	print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	print $logFile "  Generating publication tables for all peptides\n";
	print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	$publication -> generatePublicationTables(\%params, "id_all_pep.txt", $rawScanZeroTxt, $publicationDir);
	print "  Generating publication tables for unique proteins\n";
	print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	print $logFile "  Generating publication tables for unique proteins\n";
	print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	$publication -> generatePublicationTables(\%params, "id_uni_prot.txt", $rawScanZeroTxt, $publicationDir);
	print "  Generating publication tables for all proteins\n";
	print "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	print $logFile "  Generating publication tables for all proteins\n";
	print $logFile "    It may take several minutes when there are lots of PSMs, peptides and proteins\n";
	$publication -> generatePublicationTables(\%params, "id_all_prot.txt", $rawScanZeroTxt, $publicationDir);
}

print "\n";
print $logFile "\n";

##################################
## Organization of output files	##
##################################
print "  Organizing output files and directories\n";
print $logFile "  Organizing output files and directories\n";

## Delete intermediate statistical testing result files
my $statFiles = $saveDir."/"."*stat.txt";
system ("rm -rf $statFiles > /dev/null 2>&1");

## Move intermediate files
my $normFiles = $saveDir."/norm*.txt";
system ("mv $normFiles $intermediateDir > /dev/null 2>&1");
foreach (keys %{$params{'comparisonNames'}}) {
	my $compName = $params{'comparisonNames'}{$_};
	my $compFiles = $saveDir."/".$compName."*.txt";
	system ("mv $compFiles $intermediateDir > /dev/null 2>&1");
}

## Move files to 'simplified publication' directory
if (basename($idTxt) eq "ID.txt") {
	system ("cp '$publicationDir/id_uni_prot_quan.xlsx' '$simplePublicationDir' > /dev/null 2>&1");
	system ("cp '$publicationDir/id_uni_prot_zero_psm_only_quan.xlsx' '$simplePublicationDir' > /dev/null 2>&1");
} elsif (basename($idTxt) eq "IDmod.txt") {
	system ("cp '$publicationDir/id_uni_pep_quan.xlsx' '$simplePublicationDir' > /dev/null 2>&1");
	system ("cp '$publicationDir/id_uni_pep_zero_psm_only_quan.xlsx' '$simplePublicationDir' > /dev/null 2>&1");
	my $combFile = $publicationDir."/id_uni_pep_prot_combined_quan.xlsx";
	if (-e $combFile) {
		system ("cp '$combFile' '$simplePublicationDir' > /dev/null 2>&1");
	}
	my $siteFile = $publicationDir."/id_uni_site_quan.xlsx";
	if (-e $siteFile) {
		system ("cp '$siteFile' '$simplePublicationDir' > /dev/null 2>&1");
	}
}

($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
print "\n  jump -q finished: ";
printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year+1900, $mon+1, $mday, $hour, $min, $sec;
print $logFile "\n  jump -q finished: ";
printf $logFile "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year+1900, $mon+1, $mday, $hour, $min, $sec;

##################
## Subroutines	##
##################
sub printCorrMatrix {
	my ($corr, $params, $logFile) = @_;
	my @reporters = split(/\;/, $$params{'tmt_reporters_used'});	 
	my $nReporters = scalar(@reporters);
	print "          \t", join("\t", @reporters), "\n";
	print $logFile "          \t", join("\t", @reporters), "\n";	
	for (my $i = 0; $i < $nReporters; $i++) {
		print "  $reporters[$i]\t";
		print $logFile "  $reporters[$i]\t";
		for (my $j = 0; $j < $nReporters; $j++) {
			printf "%.2f\t", $$corr[$i * $nReporters + $j]; 
			printf $logFile "%.2f\t", $$corr[$i * $nReporters + $j];
		}
		print "\n";
		print $logFile "\n";
	}
	print "\n";
	print $logFile "\n";
}

sub uniq {
    my %seen = ();
    grep { not $seen{$_}++ } @_;
}

sub help {
	my ($value) = @_;
	if ($value == 0){
		print "\n";
		print "     Usage: jump -q jump_q.params \n";
		print "\n";
	}
	exit;
}


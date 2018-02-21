#!/bin/env perl

use Getopt::Long;

use Cwd;
use Cwd 'abs_path';
use File::Basename;

our $version = 1.13.0;

my $Bin=$ENV{"JUMP_L_LIB"};
use lib $ENV{"JUMP_L_LIB"};
use Clone qw(clone);
use Spiders::ProcessingMzXML;
use Spiders::Params;
use List::Util qw(shuffle);
use Parallel::ForkManager;
use Spiders::LSFQueue;
use Spiders::SGEQueue;

my ($help,$parameter,$raw_file);
GetOptions('-help|h'=>\$help,
	'-p=s'=>\$parameter,
);

my $progname = $0;
# remove path from our name
$progname =~ s@(.*)/@@i;

usage() if ($help || !defined($parameter));
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); 
my $Hydrogen_mass = 1.007825032;

##################################################
## Read a parameter file and initialization	##
##################################################
my $p = Spiders::Params->new('-path'=>$parameter);
my $params = $p->parse_param();
$params->{'Mn_Mn1'} = 0.5;
$params->{'M_M1'} = 0.3;
$params->{'low_mass'} = 57;
$params->{'high_mass'} = 187;
my $IDmod = $params->{'IDmod'};
my $Outfile = $params->{'Output'};
my $dir = (split(/\//, $IDmod))[-2];
$dir =~ s/^sum/loc/;
$dir = getcwd()."/".$dir;
if (!-e $dir) {
	system(qq(mkdir $dir >/dev/null 2>&1));
}
system(qq(cp $parameter $dir));
my $JUMPLRES = $dir."/.ID.lscore";
open(JUMPLRES,">", $JUMPLRES);
my $logFile = $dir."/jump_l.log";
open(LOGFILE,">", $logFile);
print "\n\n  Start: ";
printf "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
print LOGFILE "\n\n  Start: ";
printf LOGFILE "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
print "  Initializing jump -l program\n\n";
print LOGFILE "  Initializing jump -l program\n\n";

my $queue;
if( $params->{'cluster'} == 1 && $params->{'cluster_type'} == "LSF" ) {
    $queue = Spiders::LSFQueue->new();
}
elsif( $params->{'cluster'} == 1 && $params->{'cluster_type'} == "LSF" ) {
    $queue = Spiders::SGEQueue->new();
}
else {
    die "cluster mode selected but queueing system $parars->{'cluster_type'} unknown";
}

##########################
## Read IDmod.txt file	##
##########################
my ($frac_scan,$databaseHeader) = read_IDmod($IDmod);
print "\n";
print LOGFILE "\n";

##########################################################################################
## For each fraction, extract scan information from a .mzXML file and create .dta files	##
##########################################################################################
my $MAX_PROCESSES = 32;
print "  Extracting m/z and generating dta files for each fraction \n";
print LOGFILE "  Extracting m/z and generating dta files for each fraction \n";
if ($params -> {'cluster'} eq '1') {
	my $nTotalJobs = 0;
	my %jobIDs;
	foreach my $fraction (keys %$frac_scan) {
		my $basename = basename($fraction);
		my $new_path = $dir."/$basename";
		system(qq(mkdir $new_path >/dev/null 2>&1));
		system(qq(cp $parameter "$new_path/" >/dev/null 2>&1));	
		my $jobName = "Extraction_".$nTotalJobs;
		my$job = $queue->submit_job($new_path,$jobName,"Extraction_runshell.pl -fraction $fraction -outdir $dir -IDmod $IDmod -parameter $parameter\n\n");

		$jobIDs{$job} = 1;
		$nTotalJobs++;
		print "\r  $nTotalJobs jobs are submitted";
	}
	print "\n  You submitted $nTotalJobs job(s) for data extraction \n";
	print LOGFILE "  You submitted $nTotalJobs job(s) for data extraction \n";
	CheckJobStat($nTotalJobs, \%jobIDs, $queue);
} elsif ($params -> {'cluster'} eq '0') {
	my $nTotalJobs = 0;
	my $pm = new Parallel::ForkManager($MAX_PROCESSES);
	foreach my $fraction (keys %$frac_scan) {	
		$nTotalJobs++;
		print "\r  $nTotalJobs job(s) submitted";
		$pm -> start and next;
		my $basename = basename($fraction);
		my $new_path = $dir."/$basename";
		system(qq(mkdir $new_path >/dev/null 2>&1));
		system(qq(cp $parameter "$new_path/" >/dev/null 2>&1));
		my $jobName = "Extraction_" . $nTotalJobs;
		open (JOB, ">", "$new_path/$jobName.sh") or die "Cannot creat a job file\n";
		print JOB "#!/bin/bash\n";
		print JOB "#\$ -N $jobName\n";
		print JOB "#\$ -e $new_path/$jobName.e\n";
		print JOB "#\$ -o $new_path/$jobName.o\n";
		print JOB "Extraction_runshell.pl -fraction $fraction -outdir $dir -IDmod $IDmod -parameter $parameter\n\n";
		close (JOB);
		system(qq(sh "$new_path/$jobName.sh" > /dev/null 2>&1));
		$pm -> finish;
	}
	$pm -> wait_all_children;
	print "\n  $nTotalJobs job(s) finished\n";
	print LOGFILE "\n  $nTotalJobs job(s) finished\n";
}

##########################################################################################
## For each peptide, calculate localization scores of all possible modified peptides	##
##########################################################################################
print "\n\n  Calculating localization scores\n";
print LOGFILE "\n\n  Calculating localization scores";
if ($params -> {'cluster'} eq '1') {
	my $nEntriesPerJob = 100;
	my $nTotalJobs = 0;
	my $maxJobs = 100;
	my %jobIDs;
	foreach my $frac (keys %$frac_scan) {
		my $basename = basename($frac);
		my $new_path = $dir."/$basename";
		my @outfiles = keys %{$frac_scan->{$frac}};
		my $nEntries = scalar(@outfiles);
		my $nJobs = int($nEntries / $nEntriesPerJob) + 1;
		if ($nJobs > $maxJobs) {
			$nEntriesPerJob = int($nEntries / $maxJobs) + 1;
			$nJobs = int($nEntries / $nEntriesPerJob) + 1;
		}      
		for (my $i = 0; $i < $nJobs; $i++) {
			my $jobName = "Job_PTM_".$nTotalJobs;
			my $cmd = "";
			for (my $j = 0; $j < $nEntriesPerJob; $j++) {
				my $k = $nEntriesPerJob * $i + $j;
				last if ($k >= $nEntries);
				my $queryOutfile =  $outfiles[$k];
				my $queryPeptide = "\"" . $frac_scan->{$frac}->{$queryOutfile}->{'peptide'} . "\"";
				$cmd .= "JUMPl_runshell.pl -fraction $frac -outdir $new_path -scan $queryOutfile -peptide $queryPeptide -parameter $parameter\n\n";
			}

		        my $job = $queue->submit_job($new_path,$jobName,$cmd);
			$jobIDs{$job} = 1;
			$nTotalJobs++;
			print "\r  $nTotalJobs jobs are submitted";
		}
	}
	print "\n  You submitted $nTotalJobs job(s) for local scoring \n";
	print LOGFILE "\n  You submitted $nTotalJobs job(s) for local scoring \n";
	CheckJobStat($nTotalJobs, \%jobIDs, $queue);
} elsif ($params -> {'cluster'} eq '0') {
	my $nEntriesPerJob = 100;
	my $nTotalJobs = 0;
	my $maxJobs = 100;
	$pm = new Parallel::ForkManager($MAX_PROCESSES);
	foreach my $frac (keys %$frac_scan) {
		my $basename = basename($frac);
		my $new_path = $dir."/$basename";
		my @outfiles = keys %{$frac_scan->{$frac}};
		my $nEntries = scalar(@outfiles);
		my $nJobs = int($nEntries / $nEntriesPerJob) + 1;
		if ($nJobs > $maxJobs) {
			$nEntriesPerJob = int($nEntries / $maxJobs) + 1;
			$nJobs = int($nEntries / $nEntriesPerJob) + 1;
		}
		for (my $i = 0; $i < $nJobs; $i++) {
			$nTotalJobs++;
			print "\r  $nTotalJobs job(s) submitted (some of them may be finished)";
			my $jobName = "Job_PTM_".$nTotalJobs;
			$pm -> start and next;
			open (JOB, ">", "$new_path/$jobName.sh") or die "Cannot creat a job file\n";
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N $jobName\n";
			print JOB "#\$ -e $new_path/$jobName.e\n";
			print JOB "#\$ -o $new_path/$jobName.o\n";
			for (my $j = 0; $j < $nEntriesPerJob; $j++) {
				my $k = $nEntriesPerJob * $i + $j;
				last if ($k >= $nEntries);
				my $queryOutfile =  $outfiles[$k];
				my $queryPeptide = "\"" . $frac_scan->{$frac}->{$queryOutfile}->{'peptide'} . "\"";
				print JOB "JUMPl_runshell.pl -fraction $frac -outdir $new_path -scan $queryOutfile -peptide $queryPeptide -parameter $parameter\n\n";
			}
			close (JOB);
			system(qq(sh "$new_path/$jobName.sh" > /dev/null 2>&1));
			$pm -> finish;
		}		
	}
	$pm -> wait_all_children;
	print "\n  $nTotalJobs job(s) finished\n";
	print LOGFILE "\n  $nTotalJobs job(s) finished\n";
}

sub CheckJobStat {
	my ($nJobs, $jobIDs,$queue) = @_;
	my $nFinishedJobs = 0;
	my $jobInfo = 1;
	my ($username) = getpwuid($<);
	$| = 1;
	while($jobInfo) {
		my @jobStatusArray = @{$queue->get_running_jobs($jobIDs)};#split(/\n/,$jobStatus);
		if (@jobStatusArray) {	# i.e. There are some jobs in the queue (may include other jobs like searching)
			my @jobIDsInQueue;
			foreach (@jobStatusArray) {
				my ($jobID) = ($_ =~ /([0-9]+)\s*/);
				if (defined $$jobIDs{$jobID}) {
					push (@jobIDsInQueue, $jobID);
				}
			}
			if (@jobIDsInQueue) { # i.e. There are jobs of interest in the queue
				$nFinishedJobs = $nJobs - scalar(@jobIDsInQueue);
				print "\r  $nFinishedJobs jobs are finished";
				if ($nFinishedJobs == $nJobs) {
					$jobInfo = 0;
				} else {
					$jobInfo = 1;
				}
			} else {	# i.e. There are jobs in the queue, but all jobs of interest are finished
				print "\r  $nJobs jobs are finished";
				$jobInfo = 0;
			}
		} else {	# i.e. There's no job in the queue
			print "\r  $nJobs jobs are finished";
			$jobInfo = 0;
		}
		sleep(5);
	}
	$| = 0;
}

##########################################################
## Create .ID.lscore file (a hidden intermediate file)	##
##########################################################
print "\n";
print JUMPLRES "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;pos;precursor_peak_intensity_percentage;JUMPl_peptide;JUMPl_site;JUMPl_score\n";
foreach my $fraction (keys %$frac_scan) {
	my $basename = basename($fraction);
	my $new_path = $dir."/$basename";
	foreach my $scan (keys %{$frac_scan->{$fraction}}) {
		print "\r  Summarizing scans: $scan";
		open(INFILE,"$new_path/$scan.out");
		my $data = <INFILE>;
		chomp $data;
		my @data_array = split(/\t/,$data);

		foreach my $line (@{$frac_scan->{$fraction}->{$scan}->{'ALL'}}) {
			my ($AAstartPosition) = (split(/\;/, $line))[14] =~ /^AA(\d+)to/;
			my $Jumpl_result = join(';',@data_array[3..$#data_array]);
			print JUMPLRES $line,";",$Jumpl_result,"\n";
		}
	}
}
print "\n\n";
print LOGFILE "\n\n";
close(JUMPLRES);

##########################################################
## Retrieve the information from .ID.lscore file and	##
## create hashes containing modification site scores	##
##########################################################
my $inputFile = $JUMPLRES;
my $modSymbols = "#%*";
my $tol = $params->{'pertide_score_tolerance'};
open (IDmod, "<", $inputFile) or die "Cannot open $inputFile\n";
my $header = <IDmod>;
my %pepScores;
my %protScores;

print "  Retrieving information from the summarized file (may take a while)\n";
print LOGFILE "  Retrieving information from the summarized file (may take a while)\n";
my $nLines = 0;
while (<IDmod>) {
	$nLines++;
	print "\r  Gathering information from $nLines PSMs";
	chomp;
	my $line = $_;
	my @elems = split(/;/, $line);
	my $peptide = $elems[0];
	my $protein = $elems[1];
	my $outfile = $elems[2];
	
	## NOTE: Assumed that modification site information starts from the 17th column of a table
	## Check function for the existence of modification site information
	if (!defined $elems[16]) {
		print "  No modification site information\n";
		print "  Check the following entry\n";
		print "  $line\n";
		print LOGFILE "  No modification site information\n";
		print LOGFILE "  Check the following entry\n";
		print LOGFILE "  $line\n";
		exit;
	}

	## Generate hashes containing site scores at the peptide and protein levels
	$pepScores{$peptide}{$outfile}{'info'} = $line;
	for (my $i = 16; $i < scalar(@elems); $i = $i + 2) {
		## Find the site location in the protein
		my ($AAStartPos) = ($elems[14] =~ /AA(\d+)to/);
		my ($modAA, $pepPos) = ($elems[$i] =~ /([A-Z])(\d+)/);
		my $protSite = $modAA.($AAStartPos + $pepPos - 1);
		## Site score(s) at the peptide level
		$pepScores{$peptide}{$outfile}{'sites'}{$elems[$i]} = $elems[$i + 1];
		$pepScores{$peptide}{$outfile}{'proteins'}{$protein}{$elems[$i]} = $protSite; 
		## Site score(s) at the protein level
		if (!defined $protScores{$protein}{$protSite}) {
			$protScores{$protein}{$protSite} = $elems[$i + 1];
		} else {
			if ($elems[$i + 1] > $protScores{$protein}{$protSite}) {
				# Update the site score at the protein level
				# Take the highest site score value from the corresponding peptide scores
				$protScores{$protein}{$protSite} = $elems[$i + 1];
			}
		}
	}
}
close (IDmod);
print "\n";
print LOGFILE "  Gathered information from $nLines PSMs\n";

##########################################################
## Update the modification site(s) based on site scores	##
##########################################################
my $nRandom = 0;
my $nTotPeptides = scalar(keys %pepScores);
my $nPeptides = 0;
foreach my $peptide (keys %pepScores) {
	$nPeptides++;
	print "\r  Updating the modification-site scores in $nPeptides peptides (out of $nTotPeptides peptides)";

	## 1. Count the number and position(s) of modification site(s) in the original peptide
	my $origPeptide = $peptide;
	$origPeptide = (split(/\./, $origPeptide))[1];
	$origPeptide =~ s/\@//g;	# Remove methionine modification symbol (@)
	my $nOrigSites = 0;
	my @origSites;
	while ($origPeptide =~ /[$modSymbols]/g) {
		## n : the number of modification sites in the original peptide (=$nSites)
		$nOrigSites++;
		my $pos = pos ($origPeptide) - $nOrigSites;
		my $AA = substr($origPeptide, $pos + $nOrigSites - 2, 1);
		push (@origSites, $AA.$pos);
	}

	## 2. Go into each PSM and update the modification site(s)
	while (my ($outfile, $outfileHash) = each %{$pepScores{$peptide}}) {
		my $nSites = $nOrigSites;
		my $nSTYs = scalar(keys %{$$outfileHash{'sites'}});
		## When there's no possible modification sites other than the identified one(s),
		## we can use the original site(s)/score(s) as a localization result
		## e.g. peptide = ABCDES#FGHI
		if ($nSTYs == $nSites) {
			$$outfileHash{'selectedSites'} = \%{$$outfileHash{'sites'}};
			$$outfileHash{'selectedLevel'} = "peptide";
			next;
		}
		## When there are possible modification sites other than the identified one(s)
		## e.g. peptide = ABCDS#EFGSHITJKL
		my @pepSiteScores;
		foreach my $site (keys %{$$outfileHash{'sites'}}) {
			push (@pepSiteScores, $$outfileHash{'sites'}{$site});
		}
		@pepSiteScores = sort {$b <=> $a} @pepSiteScores;

		##########################################################################
		## 1. Choose sites as many as possible based on peptide-level scores	##
		##########################################################################
		## Select possible site-scores from @pepSiteScores
		for (my $i = $nSites - 1; $i < scalar(@pepSiteScores) - 1; $i++) {
			if ($pepSiteScores[$i] > $pepSiteScores[$i + 1] + $tol) {
				@pepSiteScores = @pepSiteScores[0..$i];
				last;
			}
		}
		## Clearly, all sites can be clearly selected based on peptide-level scores
		my $nSelected = 0;
		if (scalar(@pepSiteScores) == $nSites) {
			for (my $i = 0; $i < $nSites; $i++) {
				foreach my $site (keys %{$$outfileHash{'sites'}}) {
					if ($$outfileHash{'sites'}{$site} == $pepSiteScores[$i]) {
						$$outfileHash{'selectedSites'}{$site} = $pepSiteScores[$i];
						$$outfileHash{'selectedLevel'} = "peptide";
						$nSelected++;
					}
				}
			}
			next;
		## Some sites are ambiguous in terms of peptide-level scores
		## Select as many sites as possible at the peptide level
		} elsif (scalar(@pepSiteScores) > $nSites) {
			for (my $i = 0; $i < $nSites; $i++) {
				if ($pepSiteScores[$i] > $pepSiteScores[$nSites] + $tol) {
					foreach my $site (keys %{$$outfileHash{'sites'}}) {
						if ($$outfileHash{'sites'}{$site} == $pepSiteScores[$i]) {
							next if (defined $$outfileHash{'selectedSites'}{$site});
							$$outfileHash{'selectedSites'}{$site} = $pepSiteScores[$i];
							$nSelected++;
						}
					}
				}
			}
		} else {
			print "  Fewer peptide-level site scores than the modification sites\n";
			print "  Check $peptide in $outfile\n";
			print LOGFILE "  Fewer peptide-level site scores than the modification sites\n";
			print LOGFILE "  Check $peptide in $outfile\n";
			exit;
		}
		if ($nSelected == $nSites) {
			## Exit if all sites are determined based on peptide-level scores
			next;
		} elsif ($nSelected > 0 && $nSelected < $nSites) {
			## If some sites are determined based on peptide-level scores, those sites should not be considered hereafter
			$nSites = $nSites - $nSelected;
			for (my $i = 0; $i < $nSelected; $i++) {
				shift @pepSiteScores;
			}
			## Go to the next step; consider protein-level site scores and choose modification sites
		} elsif ($nSelected > $nSites) {
			print "  More peptide-level site scores are selected than the modification sites\n";
			print "  Check $peptide in $outfile\n";
			print LOGFILE "  More peptide-level site scores are selected than the modification sites\n";
			print LOGFILE "  Check $peptide in $outfile\n";
			exit;
		}

		##################################################################################
		## 2. Since some sites have the same peptide-level scores, 			##
		##    introduce protein-level scores in order to choose modification sites	##
		##################################################################################
		## Find the peptide modification sites to be considered
		## In @pepSiteScores, very low scored sites are already eliminated
		## 1 modification site: sites with scores >= (highest score - tolerance) will be considered hereafter
		## 2 modification sites: sites with scores >= (second highest score - tolerance) will be considered hereafter
		## 3 modification sites: sites with scores >= (third highest score - tolerance) will be considered hereafter
		## etc.
		my @pepPossibleSites;
		foreach my $site (keys %{$$outfileHash{'sites'}}) {
			if ($$outfileHash{'sites'}{$site} >= $pepSiteScores[$nSites - 1] - $tol) {
				next if (defined $$outfileHash{'selectedSites'}{$site});
				push (@pepPossibleSites, $site);
			}
		}
		## Protein-level scores of those modification sites
		my %protSiteScores;
		foreach my $protein (keys %{$$outfileHash{'proteins'}}) {
			foreach my $pepPossibleSite (@pepPossibleSites) {
				my $protSite = $$outfileHash{'proteins'}{$protein}{$pepPossibleSite};
				if (!defined $protSiteScores{$pepPossibleSite}) {
					$protSiteScores{$pepPossibleSite} = $protScores{$protein}{$protSite};
				} else {
					if ($protSiteScores{$pepPossibleSite} < $protScores{$protein}{$protSite}) {
						$protSiteScores{$pepPossibleSite} = $protScores{$protein}{$protSite};
					}
				}
			}
		}
		my @protSiteScoresArray;
		foreach my $site (keys %protSiteScores) {
			push (@protSiteScoresArray, $protSiteScores{$site});
		}
		@protSiteScoresArray = sort {$b <=> $a} @protSiteScoresArray;
		
		## Choose sites whose protein-level scores are much greater than others
		$nSelected = 0;
		for (my $i = 0; $i < $nSites; $i++) {
			if ($protSiteScoresArray[$i] > $protSiteScoresArray[$i + 1]) {
				foreach my $site (keys %protSiteScores) {
					if ($protSiteScores{$site} == $protSiteScoresArray[$i]) {
						$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
						$nSelected++;
						## If some sites are determined based on protein-level scores, those sites should not be considered hereafter
						delete $protSiteScores{$site};
					}
				}
			}
		}
		if ($nSelected == $nSites) {
			## If all sites are determined based on protein-level scores
			## 1. Assign alternative sites
			## 2. Go to the next outfile or peptide
			foreach my $site (keys %protSiteScores) {
				$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
			}
			next;
		} elsif ($nSelected > 0 && $nSelected < $nSites) {
			## If some sites are determined based on protein-level scores, those sites should not be considered hereafter
			$nSites = $nSites - $nSelected;
			for (my $i = 0; $i < $nSelected; $i++) {
				shift @protSiteScoresArray;
			}
			## Go to the next step; look for "SP" sites and choose them if exist
		} elsif ($nSelected > $nSites) {
                        print "  More peptide- and protein-level site scores are selected than the modification sites\n";
                        print "  Check $peptide in $outfile\n";
                        print LOGFILE "  More peptide- and protein-level site scores are selected than the modification sites\n";
                        print LOGFILE "  Check $peptide in $outfile\n";
			exit;
		}

		##################################################################################
		## 3. Some sites still have the same scores even at the protein-level			##
		##    We have to apply our own rules to select reasonable modification sites	##
		##################################################################################
		## Rule 1. Choose "SP" sites
		$nSelected = 0;
		my @SPsites;
		my $checkPeptide = $origPeptide;
		$checkPeptide =~ s/[$modSymbols]//g;
		foreach my $site (keys %protSiteScores) {
			my ($pos) = $site =~ /[A-Z](\d+)/;
			if (substr($checkPeptide, $pos - 1, 2) eq "SP") {
				push (@SPsites, $site);
			}
		}
		if (@SPsites) {
			if (scalar(@SPsites) > $nSites) {
				## Check if any "SP" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $SPsite (@SPsites) {
						if ($origSite eq $SPsite) {
							$nOverlap++;
						}
					}
				}
				if ($nOverlap > $nSites) {
					## If there are more SP sites (overlapped with the original modification sites) than $nSites,
					## we need to randomly select $nSites SP sites
					## Example,
					## Original peptide = SKLS#PS#PSLR
					## Candidate sites = S1, S4 and S6
					## PepSiteScores = (65.26, 65.26 and 63.80) for (S1, S4 and S6)
					## ProtSiteScores = (81.39, 80.40, 80.40) for (S1, S4 and S6)
					## - S1 is already selected based on protein-site scores
					## - One site need to be determined among S4 and S6
					## - Both of them are "SP" sites and overlapped with the original modification sites
					## - One of them needs to be "somehow" selected
					## - Then, go to the next outfile or peptide
					@SPsites = shuffle(@SPsites);
					for (my $i = 0; $i < scalar(@SPsites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@SPsites); $i++) {
							## Accept the "SP" site(s) matched to the original modification site(s)
							if ($origSite eq $SPsites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@SPsites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@SPsites = shuffle(@SPsites);
					for (my $i = 0; $i < scalar(@SPsites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "SP" sites
				foreach my $site (@SPsites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				# Go to the next outfile if all sites have been determined
				## Otherwise, go to next step; look for and choose "S" sites
				next if ($nSites == 0);
			}
		}
		
		## Rule 2. Choose "S" site(s)
		$nSelected = 0;
		my @Ssites;
		foreach my $site (keys %protSiteScores) {
			if ($site =~ /^S/) {
				push (@Ssites, $site);
			}
		}
		if (@Ssites) {
			if (scalar(@Ssites) > $nSites) {
				## Check if any "SP" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $Ssite (@Ssites) {
						if ($origSite eq $Ssite) {
							$nOverlap++;
						}
					}
				}
				if ($nOverlap > $nSites) {
					## If there are more S sites (overlapped with the original modification sites) than $nSites,
					## we need to randomly select $nSites SP sites
					## Example,
					## Original peptide = SRS#TS#EPEEAELSLSLAR
					## Candidate sites = S1, S3, T4 and S5
					## PepSiteScores = (49.91, 49.91, 49.90, 49.90) for (S1, S3, T4 and S5)
					## ProtSiteScores = (92.03, 49.99, 49.99, 49.99) for (S1, S3, T4 and S5)
					## - S1 is already selected based on protein-site scores
					## - One site need to be determined among S3 and S5
					## - Both of them are "S" sites and overlapped with the original modification sites
					## - One of them needs to be "somehow" selected
					## - Then, go to the next outfile or peptide
					@Ssites = shuffle(@Ssites);
					for (my $i = 0; $i < scalar(@Ssites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@Ssites); $i++) {
							## Accept the "S" site(s) matched to the original modification site(s)
							if ($origSite eq $Ssites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@Ssites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@Ssites = shuffle(@Ssites);
					for (my $i = 0; $i < scalar(@Ssites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "S" sites
				foreach my $site (@Ssites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				# Go to the next outfile if all sites have been determined
				## Otherwise, go to next step; look for and choose "T" sites
				next if ($nSites == 0);
			}
		}
		
		## Rule 3. Choose "T" site(s)
		$nSelected = 0;
		my @Tsites;
		foreach my $site (keys %protSiteScores) {
			if ($site =~ /^T/) {
				push (@Tsites, $site);
			}
		}
		if (@Tsites) {
			if (scalar(@Tsites) > $nSites) {
				## Check if any "T" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $Tsite (@Tsites) {
						if ($origSite eq $Tsite) {
							$nOverlap++;
						}	
					}
				}
				if ($nOverlap > $nSites) {
					@Tsites = shuffle(@Tsites);
					for (my $i = 0; $i < scalar(@Tsites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@Tsites); $i++) {
							## Accept the "T" site(s) matched to the original modification site(s)
							if ($origSite eq $Tsites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@Tsites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@Tsites = shuffle(@Tsites);
					for (my $i = 0; $i < scalar(@Tsites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "T" sites
				foreach my $site (@Tsites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				# Go to the next outfile if all sites have been determined
				## Otherwise, go to next step; look for and choose "T" sites
				next if ($nSites == 0);
			}
		}
		
		## Rule 4. Choose "Y" site(s)
		$nSelected = 0;
		my @Ysites;
		foreach my $site (keys %protSiteScores) {
			if ($site =~ /^Y/) {
				push (@Ysites, $site);
			}
		}
		if (@Ysites) {
			if (scalar(@Ysites) > $nSites) {
				## Check if any "Y" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $Ysite (@Ysites) {
						if ($origSite eq $Ysite) {
							$nOverlap++;
						}
					}
				}
				if ($nOverlap > $nSites) {
					@Ysites = shuffle(@Ysites);
					for (my $i = 0; $i < scalar(@Ysites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@Ysites); $i++) {
							## Accept the "Y" site(s) matched to the original modification site(s)
							if ($origSite eq $Ysites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@Ysites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@Ysites = shuffle(@Ysites);
					for (my $i = 0; $i < scalar(@Ysites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "Y" sites
				foreach my $site (@Ysites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				if ($nSites > 0) {
					print "  Some sites still need to be determined in $peptide of $outfile\n";
					print LOGFILE "  Some sites still need to be determined in $peptide of $outfile\n";
					exit;
				}
			}
		}
	}
}
print "\n";
print LOGFILE "  Updated the modification-site scores in $nPeptides peptides (out of $nTotPeptides peptides)\n";

##########################################################
## Generate an output file (default = ID.lscore)	##
##########################################################
my $outputFile = $dir."/$Outfile";
open (OUT, ">", $outputFile);
open (IDmod, "<", $inputFile) or die "Cannot open $inputFile\n";
$header = <IDmod>;
print OUT $databaseHeader,"\n";
print OUT $header;

my $nOutputLines = 0;
while (<IDmod>) {
	$nOutputLines++;
	print "\r  Writing the result of $nOutputLines PSMs";
	chomp;
	my $line = $_;
	my @elems = split(/;/, $line);
	my $peptide = $elems[0];
	my $outfile = $elems[2];
	my ($AAstartPos) = $elems[14] =~ /AA(\d+)to/;
	
	my @origSiteInfoArray;
	for (my $i = 16; $i < scalar(@elems); $i = $i + 2) {
		my $origSiteInfo = $elems[$i].":".$elems[$i + 1];
		push (@origSiteInfoArray, $origSiteInfo);
	}
	my @siteInfoArray;
	my @sitePosArray;
	foreach my $site (keys %{$pepScores{$peptide}{$outfile}{'selectedSites'}}) {
		my ($siteAA, $sitePos) = $site =~ /([A-Z])(\d+)/; 
		my $siteScore = $pepScores{$peptide}{$outfile}{'selectedSites'}{$site};
		my $siteInfo = $siteAA.($sitePos + $AAstartPos - 1).":".$siteScore;
		push (@sitePosArray, $sitePos);
		push (@siteInfoArray, $siteInfo);
	}
	my $modPeptide = $peptide;
	$modPeptide =~ s/[$modSymbols]//g;
	my $nMethionines = 0;
	while ($modPeptide =~ /\@/g) {
		## n : the number of modification sites in the original peptide (=$nSites)
		$nMethionines++;
		my $pos = pos ($modPeptide) - $nMethionines - 2;
		push (@sitePosArray, $pos);
	}
	$modPeptide =~ s/\@//g;
	@sitePosArray = sort {$a <=> $b} @sitePosArray;
	my $nSites = 0;	
	for (my $i = 0; $i < scalar(@sitePosArray); $i++) {
		my $AA = substr($modPeptide, $sitePosArray[$i] - 1 + $nSites + 2, 1);
		if ($AA eq "S") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "#");
		} elsif ($AA eq "T") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "%");
		} elsif ($AA eq "Y") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "*");
		} elsif ($AA eq "M") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "@");
		}
		$nSites++;
	}
	my $outputLine = join(";", @elems[0..15]).";".join(",", @origSiteInfoArray).";".$modPeptide.";".join(",", @siteInfoArray);
	print OUT "$outputLine\n";
}
close (IDmod);
close (OUT);
print LOGFILE "  Wrote the result of $nOutputLines PSMs";

##################################
## Generate publication tables	##
##################################
print "\n\n  Generating a publication table\n";
print LOGFILE "\n\n  Generating a publication table\n";
my $modSites = "STY";
# my $IDmod='IDmod_updated.txt';
my $IDmod = $dir."/IDmod.txt";
generate_pub_tables($outputFile,$modSites,$IDmod);

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
print "\n  Finished jump -l running\n";
printf "  %4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
print LOGFILE "\n  Finished jump -l running\n";
printf LOGFILE "  %4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;




#---------------------------------------------------------------------------
sub generate_pub_tables {
	my $scoreMinCut=30;
	my $scoreMedCut=60;
	my $scoreHighCut=90;
	my ($IDmod,$modSites,$output)=@_;

	my (%pephash,%prohash,$line,%out2pep,%headers,%pep2pro);
	open(IN,"$IDmod");
	$headers{IDmod}[1]=<IN>;
	$headers{IDmod}[2]=<IN>;

	# build %pephash{jumpf/jumpl}
	while(<IN>) {
		chomp;
		my @t=split(/\;/,$_);
		my ($outfile,$pep,$pro,$pos,$pepl,@AscoreSite,@Ascore,$ptmSites,$uniPep);
		($outfile,$pep,$pro,$pos,$pepl,$ptmSites,$uniPep)=($t[2],$t[0],$t[1],$t[14],$t[$#t-1],$t[$#t],$t[12]);
		next if ($pro =~ m/Decoy/);

		# mimic old format: assign @AscoreSite & @Ascore
		# S253:97.70,S248:100.00
		my @lsites=split /\,/,$ptmSites;
		my $k=0;
		foreach my $tmp (@lsites) {
			$k++;
			($AscoreSite[$k],$Ascore[$k])=split /\:/,$tmp;
		}

		# mimic old format: for $pepl, rm flanking AAs
		chop($pepl); chop($pepl); $pepl=reverse($pepl);
		chop($pepl); chop($pepl); $pepl=reverse($pepl);

		# %out2pep
		$out2pep{$outfile}{lines}{$_}='';

		# %pephash
		$pep =~ s/\.//g; 
		my $rightAA=chop($pep); $pep=reverse($pep); 
		my $leftAA=chop($pep); $pep=reverse($pep);

		$pephash{jumpf}{$pep}{'protein'}{$pro}{$pos}='';
		$pephash{jumpl}{$pepl}{'protein'}{$pro}{$pos}='';

		$pep2pro{$pepl}{$pro}=$uniPep;

		# %out2pep
		$out2pep{$outfile}{jumpf}=$pep;
		$out2pep{$outfile}{jumpl}=$pepl;
		$out2pep{$outfile}{rightAA}=$rightAA;
		$out2pep{$outfile}{leftAA}=$leftAA;

		my @mod=split(//,$modSites);
		my %mods; foreach my $s (@mod) { $mods{$s}=''; }

		# mod site pos for jump -f pep
		if (!defined($pephash{jumpf}{$pep}{modsite})) {
			my $p=$pep;
			while ($p =~ /[\%\@\#\*\^\~\$]/) {
				my $s=substr($p,$-[0]-1,1);
				if (defined($mods{$s})) { 
					my $pos=$-[0]-1; 
					$pephash{jumpf}{$pep}{'modsite'}{$pos}{score}=0; 
					$pephash{jumpf}{$pep}{'modsite'}{$pos}{AA}=$s; 
				}
				$p =~ s/[\%\@\#\*\^\~\$]//;
			}
		}

		# mod site pos for jump -l pep
		if (!defined($pephash{jumpl}{$pepl}{modsite})) {
			$pos =~ /AA(\d+)to/;
			my $startPos=$1;
			for (my $i=1; $i<=3; $i++) {
				if ( !defined($AscoreSite[$i]) or $AscoreSite[$i] eq '' ) {
					last;
				}
				$AscoreSite[$i] =~ m/(\w)(\d+)/;
				my ($modAA,$modPos)=($1,$2-$startPos);
				next unless (defined($mods{$modAA}));
				$pephash{jumpl}{$pepl}{'modsite'}{$modPos}{score}=$Ascore[$i];
				$pephash{jumpl}{$pepl}{'modsite'}{$modPos}{AA}=$modAA;
			}
		}
	}
	close IN;

	# build %prohash for jump .l peptides: only use 1st position of the 1st protein
	buildProhash($pephash{jumpl},\%prohash);

	# update %prohash for mod sites
	foreach my $outfile (keys %out2pep) {
		my $pepl=$out2pep{$outfile}{jumpl};
		my $pepf=$out2pep{$outfile}{jumpf};
		$out2pep{$outfile}{select}='l'; 
		if (!defined($pephash{select}{$pepl})) {
			$pephash{select}{$pepl}=clone($pephash{jumpl}{$pepl});
		}
		next;
	}

	# re-build %prohash  based on $pephash{select}
	undef %prohash;
	buildProhash($pephash{select},\%prohash);

	# print IDmod.txt
	open(OUT,">$output");
	print OUT "$headers{IDmod}[1]";
	if (1) {
		my @t=split /\;/,$headers{IDmod}[2];
		for (my $i=0;$i<15;$i++) {print OUT "$t[$i];";}
		print OUT "$t[15]\n";
	}
	foreach my $outfile (keys %out2pep) {
		foreach my $line (keys %{$out2pep{$outfile}{lines}}) {
			my @t=split /\;/, $line;
			if (defined($out2pep{$outfile}{select})) {
				$t[0]=($out2pep{$outfile}{select} eq 'f')?$out2pep{$outfile}{jumpf}:$out2pep{$outfile}{jumpl};
				$t[0]=join("","$out2pep{$outfile}{leftAA}\.",$t[0],"\.$out2pep{$outfile}{rightAA}");
				for (my $i=0; $i<15; $i++) { print OUT "$t[$i];"; }
				print OUT "$t[15]\n";
			}
			else { die "$outfile select not defined\n"; }
		}
	}
	close OUT;

	# calculate phospho sites
	my %counts;
	if (1) {
		my @mod=split(//,$modSites);
		foreach my $aa (@mod) {
			$counts{$aa}{low}=$counts{$aa}{med1}=$counts{$aa}{med2}=$counts{$aa}{high}=0;
		}
		my $aa='total';
		$counts{$aa}{low}=$counts{$aa}{med1}=$counts{$aa}{med2}=$counts{$aa}{high}=0;
	}
	foreach my $pro (keys %prohash) {
		foreach my $ms (keys %{$prohash{$pro}}) {
			my ($aa,$score)=($prohash{$pro}{$ms}{AA},$prohash{$pro}{$ms}{score});
			if ( $score<$scoreMinCut ) {$counts{$aa}{low}++;$counts{total}{low}++;}
			elsif ( $scoreMinCut<=$score and $score<$scoreMedCut ) 
			{$counts{$aa}{med1}++;$counts{total}{med1}++;}
			elsif ( $scoreMedCut<=$score and $score<$scoreHighCut ) 
			{$counts{$aa}{med2}++;$counts{total}{med2}++;}
			else {$counts{$aa}{high}++;$counts{total}{high}++;}
		}
	}

	print "  Phospho pep = ",scalar(keys %{$pephash{select}}),"\n";
	print "  Phospho sites:\n";
	print "  modsites\tscore<$scoreMinCut\t$scoreMinCut-$scoreMedCut\t$scoreMedCut-$scoreHighCut\t>=$scoreHighCut\tsum\n";
	print LOGFILE "  Phospho pep = ",scalar(keys %{$pephash{select}}),"\n";
	print LOGFILE "  Phospho sites:\n";
	print LOGFILE "  modsites\tscore<$scoreMinCut\t$scoreMinCut-$scoreMedCut\t$scoreMedCut-$scoreHighCut\t>=$scoreHighCut\tsum\n";
	foreach my $aa (keys %counts) {
		next if ($aa eq 'total');
		print "  $aa\t\t$counts{$aa}{low}\t$counts{$aa}{med1}\t$counts{$aa}{med2}\t$counts{$aa}{high}\t",$counts{$aa}{low}+$counts{$aa}{med1}+$counts{$aa}{med2}+$counts{$aa}{high},"\n";
		print LOGFILE "  $aa\t\t$counts{$aa}{low}\t$counts{$aa}{med1}\t$counts{$aa}{med2}\t$counts{$aa}{high}\t",$counts{$aa}{low}+$counts{$aa}{med1}+$counts{$aa}{med2}+$counts{$aa}{high},"\n";
	}
	my $aa='total';
	print "  $aa\t\t$counts{$aa}{low}\t$counts{$aa}{med1}\t$counts{$aa}{med2}\t$counts{$aa}{high}\t",$counts{$aa}{low}+$counts{$aa}{med1}+$counts{$aa}{med2}+$counts{$aa}{high},"\n";
	print LOGFILE "  $aa\t\t$counts{$aa}{low}\t$counts{$aa}{med1}\t$counts{$aa}{med2}\t$counts{$aa}{high}\t",$counts{$aa}{low}+$counts{$aa}{med1}+$counts{$aa}{med2}+$counts{$aa}{high},"\n";

	# print publication tables:
	# Get jump -f publication table paths
	my $path=$headers{IDmod}[1];
	chomp($path);
	$path =~ s/^Database=//;
	#$path =~ s/idsum_mod.db$/publications\//;
	$path =~ s/(idsum_mod.db|sum_mod.db)$/publications\//;
	my $allPep=$path; my $uniPep=$path;
	$uniPep.='id_uni_pep.txt';
	$allPep.='id_all_pep.txt';
	if (!(-e $uniPep)) { die "$uniPep not exist!!!\n"; }
	if (!(-e $allPep)) { die "$allPep not exist!!!\n"; }

	my (%proInfor,%pepInfor);

	# build %proInfor{$pro}{group/accession/description/GN/attributes(unique/all)}:
	open(IN,$allPep) || die "cannot open $allPep\n";
	$headers{allPep}[1]=<IN>;
	$headers{allPep}[2]=<IN>;
	$headers{allPep}[3]=<IN>;
	$headers{allPep}[4]=<IN>;
	$headers{allPep}[5]=<IN>;
	while(<IN>) {
		chomp;
		my ($pep0,$ProteinGroup,$pro,$Description,$GN,$PSM,$Run,$Scan,$mz,$z,$ppm,$Jscore,$dJn,$Qvalue)=split /\t/,$_;
		#Peptides        Protein Group#  Protein Accession #     Protein Description     GN      PSM#    Run#    Scan#   m/z     z       ppm     Jscore  dJn     Q-value(%)
		#R.SAAVETAS#LSPSLVPAR.Q  SJPG105434.001  sp|Q2KHT3|CL16A_HUMAN   Protein CLEC16A OS=Homo sapiens GN=CLEC16A PE=2 SV=2    CLEC16A 1       p1f17   84765   982.51044       2  -0.769628420271329       32.40   0.547549496916585       0.229417206290472
		$proInfor{$pro}{group}=$ProteinGroup;
		$proInfor{$pro}{description}=$Description;
		$proInfor{$pro}{GN}=$GN;
		$proInfor{$pro}{attributes}='all';
		$proInfor{$pro}{uniPepN}=0;
		$proInfor{$pro}{pepN}=0;

		my $pep=$pep0;
		chop($pep); chop($pep); $pep=reverse($pep);
		chop($pep); chop($pep); $pep=reverse($pep);
		if (defined($pephash{jumpf}{$pep})) { 
			$pephash{jumpf}{$pep}{Qvalue}=$Qvalue; 
		}
		#else { die "$pep in $allPep not found in IDascore.txt!!!\n"; }
	}
	close IN;

	open(IN,$uniPep) || die "cannot open $uniPep\n";
	$headers{uniPep}[1]=<IN>;
	$headers{uniPep}[2]=<IN>;
	$headers{uniPep}[3]=<IN>;
	$headers{uniPep}[4]=<IN>;
	$headers{uniPep}[5]=<IN>;

	while(<IN>) {
		chomp;
		my ($pep0,$ProteinGroup,$pro,$Description,$GN,$PSM,$Run,$Scan,$mz,$z,$ppm,$Jscore,$dJn,$Qvalue)=split /\t/,$_;
		$proInfor{$pro}{attributes}='uni';
	}
	close IN;

	# add peptide counts (both unique and total) to %proInfor{$pro}{uniPepN} & %proInfor{$pro}{pepN}
	foreach my $pep (keys  %pep2pro) {
		foreach my $pro (keys %{$pep2pro{$pep}}) {
			$proInfor{$pro}{pepN}++;
			if ( $pep2pro{$pep}{$pro} ) { $proInfor{$pro}{uniPepN}++; }
		}
	}

	# build %pepInfor
	my (%run2pep);
	foreach my $out (keys %out2pep) {
		my $pep=($out2pep{$out}{select} eq 'f')?$out2pep{$out}{jumpf}:$out2pep{$out}{jumpl};
		my $line; foreach my $l (keys %{$out2pep{$out}{lines}}) { $line=$l; last; }
		my ($Peptide,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos)=split /[\,\;]/,$line;
		next if ($Protein =~ m/Decoy/); # could cause peptide missing if T/D shared

		my @t=split /\//,$out;
		my ($run,$scan,$ppiN,$z)=split /\./,$t[$#t];
		#print "$run,$scan,$z,",int($measuredMH+0.5),"\n" if ($t[$#t] eq 'p15.112637.3.4.spout');
		$run2pep{$run}{$scan}{$z}{int($measuredMH+0.5)}=$pep;

		if (!defined($pepInfor{$pep}) or $pepInfor{$pep}{score}<$XCorr) {
			$pepInfor{$pep}{score}=$XCorr;
			$pepInfor{$pep}{mz}=($measuredMH-$Hydrogen_mass)/$z+$Hydrogen_mass;
			$pepInfor{$pep}{z}=$z;
			$pepInfor{$pep}{ppm}=$ppm;
			$pepInfor{$pep}{dCn}=$dCn;
			$pepInfor{$pep}{run}=$run;
			$pepInfor{$pep}{scan}=$scan;

			# Q-value
			my $pepf=$out2pep{$out}{jumpf};
			if (defined($pephash{jumpf}{$pepf}{Qvalue})) { $pepInfor{$pep}{Qvalue}=$pephash{jumpf}{$pepf}{Qvalue}; }
			else { die "$pepf in IDascore.txt not defined in $allPep\n"; }

			# leftAA/rightAA
			$pepInfor{$pep}{leftAA}=$out2pep{$out}{leftAA};
			$pepInfor{$pep}{rightAA}=$out2pep{$out}{rightAA};

			# protein names
			foreach $line (keys %{$out2pep{$out}{lines}})
			{
				my @t=split /[\,\;]/,$line;
				my $pro=$t[1];
				$pepInfor{$pep}{proteins}{$pro}='';
			}

			# select unique protein
			my $minGroup=1000000000;
			foreach my $pro (keys %{$pepInfor{$pep}{proteins}})
			{
				my $group=$proInfor{$pro}{group};
				if ($group =~ m/^SJPG(.*?)$/) { 
					$group=$1;
				} else {
					next;
#				 die "group pattern unmatched: $group\n"; }
				}
				if ( $minGroup>$group ) 
				{ 
					$minGroup=$group; 
					$pepInfor{$pep}{uniqPro}=$pro;
				}
			}

			# mod site
			if (defined($pephash{select}{$pep}{'modsite'}))
			{
				my $ms=min(keys %{$pephash{select}{$pep}{'modsite'}});
				$pepInfor{$pep}{lscore}=$pephash{select}{$pep}{'modsite'}{$ms}{score};

				my $pro=$pepInfor{$pep}{uniqPro};
				$pepInfor{$pep}{proGroup}=$proInfor{$pro}{group};
				my @t=keys %{$pephash{select}{$pep}{protein}{$pro}};
				my $pos=shift @t;
				if (!defined($pos)) { die "pos not defined: $pep,$pro,$pos\n"; }
				$pos =~ m/AA(\d+)toAA(\d+)/;
				my $startpos=$1;
				$ms+=$startpos; 
				$pepInfor{$pep}{ms}=$ms;
			}
			else { die "$pep not defined in pephash!!!\n"; }
		}
		$pepInfor{$pep}{PSM}++;
	}

	# print publication tables
	my $pubDir = (split(/\//, $output))[-2];
	$pubDir = $pubDir."/publications";
	unless (-e $pubDir) {
		system("mkdir $pubDir");
	}
	printPublicationTables(\%pepInfor,\%headers,1,"$pubDir/id_uni_pep.txt",\%proInfor,\%counts,\%pephash);
	printPublicationTables(\%pepInfor,\%headers,0,"$pubDir/id_all_pep.txt",\%proInfor,\%counts,\%pephash);
	printProteinPublicationTables(\%proInfor,$path,1,"$pubDir/id_uni_prot.txt",\%run2pep);
	printProteinPublicationTables(\%proInfor,$path,0,"$pubDir/id_all_prot.txt",\%run2pep);
}


#-------------------------------------------------------------------------------
sub printProteinPublicationTables
{
	my ($proInfor,$path,$unique,$output,$run2pep)=@_;
	my $input=$path;
	if ($unique) {$input.='id_uni_prot.txt';}
	else {$input.='id_all_prot.txt';}

	open(IN,$input) || die "Cannot open original protein publication table: $input!!!\n";
	#print "$input\n";
	open(OUT,">$output");
	my $header=<IN>; print OUT $header;
	$header=<IN>; print OUT $header;
	while(<IN>)
	{
		chomp;
		my @t=split /\t/,$_;
		if (defined($$proInfor{$t[1]}))
		{
			$t[5]=$$proInfor{$t[1]}{pepN};
			$t[6]=$$proInfor{$t[1]}{uniPepN};
			#if (defined($$run2pep{$t[9]}{$t[10]}{$t[12]}{int($MH+0.5)}))
			my $MH=($t[11]-$Hydrogen_mass)*$t[12]+$Hydrogen_mass;
			if (defined($$run2pep{$t[9]}{$t[10]}{$t[12]}{int($MH+0.5)}))
			{ $t[8]=$$run2pep{$t[9]}{$t[10]}{$t[12]}{int($MH+0.5)}; }
			elsif (defined($$run2pep{$t[9]}{$t[10]}{$t[12]}{int($MH+0.5)-1}))
			{ $t[8]=$$run2pep{$t[9]}{$t[10]}{$t[12]}{int($MH+0.5)-1}; }
			#else { die "$t[1] best peptide not found in run2pep hash:$t[9],$t[10],$t[11],$t[12]",int($t[11]*$t[12]+0.5),"\n"; }
		}
		else { die "protein $t[1] defined in original protein publication table but not found in proInfor hash!!!\n"; }
		for (my $i=0;$i<$#t;$i++) { print OUT "$t[$i]\t"; }
		print OUT "$t[$#t]\n";
	}
	close IN;
	close OUT;
}

sub printPublicationTables
{
	my ($pepInfor,$headers,$unique,$output,$proInfor,$counts,$pephash)=@_;
	open(OUT,">$output");
	my $totalModSites=$$counts{total}{low}+$$counts{total}{med1}+$$counts{total}{med2}+$$counts{total}{high};
	my $SModSites=$$counts{S}{low}+$$counts{S}{med1}+$$counts{S}{med2}+$$counts{S}{high};
	my $TModSites=$$counts{T}{low}+$$counts{T}{med1}+$$counts{T}{med2}+$$counts{T}{high};
	my $YModSites=$$counts{Y}{low}+$$counts{Y}{med1}+$$counts{Y}{med2}+$$counts{Y}{high};
	if ($unique)
	{
		print OUT $$headers{uniPep}[1];
		print OUT "n = ",scalar(keys %{$pepInfor})," modified peptides\n";
		print OUT "n = $totalModSites modified sites ($SModSites S, $TModSites T, $YModSites Y)\n";
		print OUT $$headers{uniPep}[4];
		chomp($$headers{uniPep}[5]);
		print OUT $$headers{uniPep}[5],"\tMod sites\tLscore\n";
	}
	else
	{
		print OUT $$headers{allPep}[1];
		print OUT "n = ",scalar(keys %{$pepInfor})," modified peptides\n";
		print OUT "n = $totalModSites modified sites ($SModSites S, $TModSites T, $YModSites Y)\n";
		print OUT $$headers{allPep}[4];
		chomp($$headers{allPep}[5]);
		print OUT $$headers{allPep}[5],"\tMod sites\tLscore\n";
	}

	for my $pep ( sort { $$pepInfor{$a}{proGroup} cmp $$pepInfor{$b}{proGroup} ||
		$$pepInfor{$a}{uniqPro} cmp $$pepInfor{$b}{uniqPro} ||
		$$pepInfor{$a}{ms} <=> $$pepInfor{$b}{ms} ||
		$$pepInfor{$b}{lscore} <=> $$pepInfor{$a}{lscore}
		} keys %{$pepInfor} )
	{
		if (!defined($$pepInfor{$pep}{uniqPro})) { die "uniqPro not define:$pep\n"; }
		for my $pro (sort { $$proInfor{$a}{group} cmp $$proInfor{$b}{group} } keys %{$$pepInfor{$pep}{proteins}})
		{
			if ($unique and $$pepInfor{$pep}{uniqPro} ne $pro) {next;}
			next if (!defined($$pephash{select}{$pep}{protein}{$pro}));

			my $mz = sprintf("%.4f", $$pepInfor{$pep}{mz});
			my $ppm = sprintf("%.2f", $$pepInfor{$pep}{ppm});
			my $score = sprintf("%.2f", $$pepInfor{$pep}{score});
			my $dCn = sprintf("%.2f", $$pepInfor{$pep}{dCn});
			my $Qvalue = sprintf("%.3f", $$pepInfor{$pep}{Qvalue});

			print OUT "$$pepInfor{$pep}{leftAA}\.$pep\.$$pepInfor{$pep}{rightAA}\t$$proInfor{$pro}{group}\t$pro\t$$proInfor{$pro}{description}\t$$proInfor{$pro}{GN}\t$$pepInfor{$pep}{PSM}\t$$pepInfor{$pep}{run}\t$$pepInfor{$pep}{scan}\t$mz\t$$pepInfor{$pep}{z}\t$ppm\t$score\t$dCn\t$Qvalue\t";

			# mod sites:
			my @t=keys %{$$pephash{select}{$pep}{protein}{$pro}};
			my $pos=shift @t;
			if (!defined($pos)) { die "pos not defined: $pep,$pro,$pos\n"; }
			$pos =~ m/AA(\d+)toAA(\d+)/;
			my $startpos=$1; 
			my @tmparray;
			for my $ms (sort {$a <=> $b} keys %{$$pephash{select}{$pep}{'modsite'}})
			{
				push @tmparray,$$pephash{select}{$pep}{'modsite'}{$ms}{AA}.($startpos+$ms);
			}
			print OUT join(",",@tmparray);
			print OUT "\t";
			undef @tmparray;
			for my $ms (sort {$a <=> $b} keys %{$$pephash{select}{$pep}{'modsite'}})
			{
				push @tmparray,$$pephash{select}{$pep}{'modsite'}{$ms}{score};
			}
			print OUT join(",",@tmparray);
			print OUT "\n";
		}
	}

	close OUT;
}

sub buildProhash
{
	my ($pephash,$prohash)=@_;
	# build %prohash for jump -l peptides: only use 1st position of the 1st protein
	# %prohash{$pro}{$modSite}{score} = max($ascore)
	# 			  {AA}    = S/T/Y
	foreach my $pep (keys %{$pephash})
	{
		my $firstPro=1;
		foreach my $pro (sort {$a cmp $b} keys %{$$pephash{$pep}{protein}})
		{
			my $firstLocation=1;
			foreach my $pos (keys %{$$pephash{$pep}{protein}{$pro}})  # peptide to protein position
			{
				$pos =~ /^AA(\d+)toAA(\d+)$/;
				my $start=$1;
				if ( $start<1 ) { print "Weired position information: $pep,$pos\n";next; }
				foreach my $ms (keys %{$$pephash{$pep}{modsite}})
				{
					my $t=$start+$ms;
					if (!defined($$prohash{$pro}{$t})) 
					{ 
						#$prohash{$pro}{$t}=$firstPro; 
						my ($aa,$score)=($$pephash{$pep}{modsite}{$ms}{AA},$$pephash{$pep}{modsite}{$ms}{score});
						$$prohash{$pro}{$t}{AA}=$aa; 
						if ( !defined($$prohash{$pro}{$t}{score}) or $$prohash{$pro}{$t}{score}<$score )
						{
							$$prohash{$pro}{$t}{score}=$score; 
						}
					}
				}
				$firstLocation=0; last; # only use first location
			}
			$$pephash{$pep}{bestprotein}=$pro;  # record best protein in %pephash
			$firstPro=0; last; # only use first protein
		}
	}
}

sub min {
    my (@vars) = @_;
    my $min;
    for (@vars) {
        $min = $_ if !$min || $_ < $min;
    }
    return $min;
}

sub string2num {
	my $st = @_;
	my @b=split(//,$st);
	my $c=0; for (@b) {$c+=ord($_);}
	return $c;
}

sub read_IDmod
{
	##########################################
	## Read ID.txt file and generate hashes ##
	##########################################
	my ($idTxt) = @_;
	
	my %frac_scan;
	
	print "  Parsing information from $idTxt\n";
	print LOGFILE "  Parsing information from $idTxt\n";

	open (INFILE, $idTxt);
	while (<INFILE>) {
		chomp;
		## Header line
		## Finally, the header includes the header of ID.txt and ASCORE header (defined above) 
		if ($_ =~ /^Peptide/) {
			$_ =~ s/\;/\,/g;
			$header = $_.",".$header;
		## Database line
		} elsif ($_ =~ /^Database/) {
			$databaseHeader = $_;
		## Peptide information line
		} else {
			my @elems = split(/\;/, $_);
			if (($elems[0] ne "") && (defined($elems[0]))) {
				#####################
				## Generate hashes ##
				#####################
				my $line = $_;
				## Key: a directory name where .out files in ID.txt file come from (used to generate and indicate ACORE.csv files)
				## Value: .out file names comes from the "key" directory
				## For example, if ID.txt file contains .out files from  three fractions, 
				## %dir2out will have three keys (directory names) representing those fractions  
				my $dir = dirname($elems[2]);
				my $file = basename($elems[2]);
				$file =~ s/(.out|.spout)//;
				push (@{$frac_scan{$dir}{$file}{'ALL'}},$line);	
				print "\r  processing: $file";
#				print LOGFILE  "\r  processing: $file";
				
				$frac_scan{$dir}{$file}{'peptide'}=$elems[0];
				$frac_scan{$dir}{$file}{'MH1'}=$elems[4];
			
				## %out2id
				## Key: .out file name (with an absolute path)
				## Value: corresponding output information in id.txt file
				## Therefore, %out2id has the keys equal to the number of unique .out files in ID.txt file
				my $idKey = $dir."/".$file;
				$_ =~ s/\;/\,/g;
				push (@{$out2id{$idKey}}, $_);
			}
		}
	}
	close INFILE;
	return (\%frac_scan,$databaseHeader);
}	

sub usage {
print <<"EOF";
Usage: $progname -p parameterfile
	or
       $progname -p parameterfile
EOF
exit 1;
}

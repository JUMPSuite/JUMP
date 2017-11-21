#!/usr/bin/perl

## General packages
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use File::Spec;
use Parallel::ForkManager;
use Time::Piece;


## Custom packages
my $libPath;
BEGIN {$libPath = getcwd;};
use lib "$libPath";
use Spiders::Params;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::Decharge;
use Spiders::BuildIndex;
use Spiders::Job;
use Spiders::Error;
use Spiders::Digestion;
use Spiders::MakeDB;
use Spiders::MassAccuracy;
use Spiders::PIP;
use Spiders::Path;
use Spiders::RankHits;
use Spiders::SpoutParser;
use Spiders::MassCorrection;

## Set the library path
my $library = $libPath;

## Get arguments and show the usage of the script if necessary
my $progname = $0;
$progname =~ s@(.*)/@@i;
my ($help, $parameter, $raw_file);
GetOptions('-help|h' => \$help,
			'-p=s' => \$parameter,);
usage() if ($help || !defined($parameter));
usage("Please verify your parameter file $parameter path is correct\n") if(!(-e $parameter));
foreach my $file (@ARGV) {
    usage("Please verify your input $file path is correct\n") if(!(-e $file));
}

unless( File::Spec->file_name_is_absolute($parameter) ) {
    $parameter = File::Spec->rel2abs($parameter);
}

print <<EOF;
################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump search using JUMP                  ****     # 
#       ****  Version 12.1.0	                      ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF

## Initialization and predefinition of some parameters
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 
my $p = Spiders::Params -> new('-path' => $parameter);
my $params = $p -> parse_param();
$params -> {'Mn_Mn1'} = 0.5;
$params -> {'M_M1'} = 0.3;
$params -> {'low_mass'} = 57;
$params -> {'high_mass'} = 187;
$params -> {'tag_coeff_evalue'} = 1;
$params -> {'pho_neutral_loss'} = 0;
database_creation();

#

## Create the path for multiple raw files
my %rawfile_hash;
print "  Using the following rawfiles:\n";
foreach my $arg (sort @ARGV) {
    my @suffixlist = ();
    push @suffixlist, ".raw";
    push @suffixlist, ".RAW";
    push @suffixlist, ".mzXML";
    if ($arg =~ /.[raw|RAW|mzXML]/) {
		print "  $arg","\n";
	}
} 

print "\t\t100 ";
print localtime->strftime('%Y%m%d %k:%M:%S');
print "\n ";


foreach my $arg (sort @ARGV) {
    my @suffixlist = ();
    push @suffixlist, ".raw";
    push @suffixlist, ".RAW";
    push @suffixlist, ".mzXML";
    if ($arg =~ /.[raw|RAW|mzXML]/) {	
		my ($filename, $directory, $suffix) = fileparse($arg, @suffixlist);	
        system(qq(mkdir $directory/$filename >/dev/null 2>&1));
        system(qq(mv $arg $directory/$filename >/dev/null 2>&1));
        my $datafile = "$directory/$filename";
        my $path = new Spiders::Path($datafile);
        my $list = $path -> make_directory_list();
        my $newdir;
        
        $newdir	= $filename . ".1";


#        if (@$list) {
#            $newdir = $path -> choose_dir_entry($list, "  Choose a .out file directory (1)", $newdir);
#        } else {
#			$newdir	= $filename . ".1";
#            $newdir = $path -> ask("  Choose a .out file directory (2)", $newdir);
#        }
        print "  Using: $newdir\n";
        $path -> add_subdir($newdir);
        my $dir =  $path -> basedir() . "/$newdir";
		my $rawfile = "$datafile/$filename";
		$rawfile_hash{$rawfile} = $dir;
		system(qq(cp -rf $parameter "$dir/jump.params" >/dev/null 2>&1));
    }
}

foreach my $raw_file (sort keys %rawfile_hash) {		
	## Get working directory
	print "  Searching data: $raw_file\n";
	my $dta_path = $rawfile_hash{$raw_file};
	if ($raw_file =~ /\.mzXML/) {
		$raw_file =~ s/\.mzXML//g;
	}

	## Conversion of .raw file(s) to .mzXML file(s)
	my $proc_raw = new Spiders::ProcessingRAW();
	$proc_raw -> set_raw_file($raw_file);
	print "  Converting .raw into .mzXML file\n";
	my $mzXML = $proc_raw->raw2mzXML();
	
	## Initialize the processing of mzXML file(s)
	print "  Extracting peaks from .mzXML\n";
	my $proc_xml = new Spiders::ProcessingMzXML();
	$proc_xml -> set_dta_path($dta_path);
	$proc_xml -> set_mzXML_file($mzXML);
                           
	## Extraction of MS1 and MS2 spectra information (assignment to hashes)
	my (%ms_hash, %msms_hash, @mz_array);
	$proc_xml -> set_parameter($params);

    printf "\n";
	print "\t\t101 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
    printf "\n\n";

    # Gathering scan information: 
    
	$proc_xml-> generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);
	my $ms1N = scalar(keys %{$ms_hash{'surveyhash'}});	## Number of MS1 scans
	my $ms2N = scalar(keys %msms_hash) - $ms1N;	## Number of MS2 scans
	printf "\n  There are %d MS and %d MS/MS in the entire run\n", $ms1N, $ms2N;
	
	printf "\n\n";
	print "\t\t102 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
    printf "\n ";

	## Mass correction of MS1 and/or MS2 spectra
	print "\n  Mass correction\n";
	my $masscorr = new Spiders::MassCorrection();
	my ($msms_hash_corrected, $mz_array_corrected) = $masscorr -> massCorrection(\%ms_hash, \%msms_hash, \@mz_array, $params);
	%msms_hash = %$msms_hash_corrected;
	@mz_array = @$mz_array_corrected;
	
	printf "\n";
	print "\t\t103 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
	printf "\n ";
	
	## Decharging/deconvolution of spectra
	print "\n  Decharging scans\n";
	my $pip = new Spiders::PIP;
	$pip -> set_parameter($params);
	$pip -> set_origmz_array(\@mz_array);
	$pip -> set_origmsms_hash(\%msms_hash);
	$pip -> set_isolation_window($params -> {'isolation_window'});
	if (!defined($params -> {'isolation_window_offset'})) {
		$params->{'isolation_window_offset'} = 0;
		print "  Warning: you did not define isolation_window_offset parameter, use default value of 0\n";
	}
	if (!defined($params -> {'isolation_window_variation'})) {
		$params->{'isolation_window_variation'} = 0;
		print "  Warning: you did not define isolation_window_variation parameter, use default value of 0\n";		
	}
	$pip -> set_isolation_offset($params -> {'isolation_window_offset'});
	$pip -> set_isolation_variation($params -> {'isolation_window_variation'});
	$pip -> set_dta_path($dta_path);
	
	printf "\n";
	print "\t\t104 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
	printf "\n ";
	
	my $PIPref = $pip -> Calculate_PIP();
	my ($charge_dist, $ppi_dist) = $pip -> changeMH_folder($PIPref);

	printf "\n";
	print "\t\t105 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
	printf "\n ";
	
	## Database search
	print "  Starting database searching\n";
	## Create jobs and distribute to nodes
	my $job = new Spiders::Job;
	$job -> set_library_path($library);
	$job -> set_dta_path("$dta_path");
	$job -> set_pip($PIPref);
	my @file_array = glob("$dta_path/*.dta");
	my $random = int(rand(100));
	if ($params -> {'second_search'} == 0) {
		$job -> create_script(0);
	} elsif ($params -> {'second_search'} == 1) {
		$job -> create_script(1);
	} else {
		print "  Please specify a right second_search parameter!!\n";
	}
		
	printf "\n";
	print "106 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
	printf "\n ";

	
	
	my $temp_file_array = runjobs(\@file_array, $dta_path, "sch_${random}");

	printf "\n\n";
	print "108 ";
	print localtime->strftime('%Y%m%d %k:%M:%S');
	printf "\n ";

	
	## Re-searching part; look for accidentally unfinished jobs and re-run them
	my $rerunN = 0;
	while (scalar(@{$temp_file_array}) > 0 && $rerunN < 3) {
		$rerunN++;
		print "\n", scalar(@$temp_file_array), " .dta files not finished! Doing re-search (rerunN = $rerunN)\n";
		my $remaining = runjobs($temp_file_array, $dta_path, "rescue_$rerunN");
		$temp_file_array = $remaining;
	}
	if (scalar(@{$temp_file_array}) > 0) {
		print "\nWarning: ", scalar(@$temp_file_array), " .dta files not finished!\n";
	}		

	## Output handling???
	
	my $Rankhit = new Spiders::RankHits();
	my $p = Spiders::Params -> new('-path'=>$parameter);
	my $params = $p -> parse_param();
	$Rankhit -> set_parameter($params);
	my $mainhash = $Rankhit -> parse_spOut_files_v5($dta_path);
	my ($sum_target, $sum_decoy, $cutoff_score) = $Rankhit -> calculate_FDR($mainhash, 0.01);
	print "\n  $sum_target targets and $sum_decoy decoys passed FDR = 1%\n";
	print "  Generating dtas file\n";
	system (qq(cat $dta_path/*.dtas > $dta_path.dtas));      # make one dtas
	system (qq(rm $dta_path/*.dtas));        # delete job-specific dtas
	print "  Generating pepXML file\n";
	my $spout = Spiders::SpoutParser -> new();
	my $outXML = $dta_path . ".pepXML"; 
	my (%seqParaHash, %outInforHash);
	$spout -> readSpsearchParam(\%seqParaHash, "$dta_path/jump.params");
	$spout -> parseSpout(\%seqParaHash, \%outInforHash, "$dta_path");
	$spout -> printPepXML(\%seqParaHash, \%outInforHash, "$dta_path", $outXML, 5);
	print "\n  Generating summary XLSX report file\n";
	$spout -> set_MS1scanNumber($ms1N);
	$spout -> set_MS2scanNumber($ms2N);
	$spout -> set_totalScan($ms1N+$ms2N);
	$spout -> set_chargeDistribution($charge_dist);
	$spout -> set_ppiDistribution($ppi_dist);
	my $report_xls = $dta_path . ".xlsx";
	my $topMS2_array = topMS2scans(\%msms_hash);	
	$spout -> printXLSX(\%seqParaHash, \%outInforHash, $dta_path, $report_xls, $topMS2_array);
	if (defined($params -> {'temp_file_removal'}) && $params -> {'temp_file_removal'} == 1) {
		delete_run_dir($dta_path);
	}	
	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 
	print "  Date: ";
	printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
}

print "  Search finished\n\n";
			
sub Create_Sort_BashFile {
	my ($params, $dta_path) = @_;
	my $FileName = "$dta_path/sort_db_" . $params->{"range"} . ".sh";
	my $cmd = join(" ","perl $dta_path/Create_Partial_Idx.pl -m", $params->{"range"}, "-job_num", $params->{"JobNum"}, "-dta_path", $dta_path, "-database_path", $dta_path,"-mass_index",$params->{"mass_index"},"-peptide_index", $params->{"peptide_index"},"-protein_index",$params->{"protein_index"},"-databasename",$params->{"databasename"},"-num_pro_per_job",$params->{"num_pro_per_job"},"-prot_job_num",$params->{"prot_job_num"},"-digestion",$params->{"digestion"},"\n");
	LaunchParallelJob($FileName, $cmd, $params -> {"GridType"}, "sort_db_" . $params -> {"range"}, $dta_path);
}

sub runjobs {
	my ($file_array, $dta_path, $job_name) = @_;
	my $MAX_PROCESSES = 32;
	my $dta_num_per_file = 10;
	my $job_num = int($#$file_array / $dta_num_per_file) + 1;

	## Set the maximum number of jobs to 4000
	if ($job_num > 4000) {
		$job_num = 4000;
		$dta_num_per_file = int($#$file_array / $job_num) + 1;
	}
	if ($params -> {'cluster'} eq '0') {
		$job_num = $MAX_PROCESSES;
		$dta_num_per_file = int($#$file_array / $job_num) + 1;
	}
	
	## Write shell scripts for parallel jobs
	for (my $i = 0; $i < $job_num; $i++) {	
		if (($i * $dta_num_per_file) > $#$file_array) {
			$job_num = $i;
			last;
		}
		open (JOB, ">", "$dta_path/${job_name}_${i}.sh") || die "can not open the job files\n";
		my $dta_file_temp = "";
		my @dta_file_arrays = ();
		my $multiple_jobs_num = 0;
		for (my $j = 0; $j < $dta_num_per_file; $j++) {
			if (($i * $dta_num_per_file + $j) <= $#$file_array) {
				$dta_file_temp .= " $$file_array[$i * $dta_num_per_file + $j]";
				push (@dta_file_arrays, $$file_array[$i * $dta_num_per_file + $j])
			}
		}
		if($params->{'Job_Management_System'} eq 'LSF') {
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q normal\n";
			print JOB "#BSUB -M 2000\n";
			print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
			print JOB "#BSUB -eo $dta_path/${job_name}_${i}.e\n";
			print JOB "#BSUB -oo $dta_path/${job_name}_${i}.o\n";
			print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $parameter -dta_path $dta_path $dta_file_temp\n";		
		} elsif($params->{'Job_Management_System'} eq 'SGE') {
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N ${job_name}_${i}\n";
			print JOB "#\$ -e $dta_path/${job_name}_${i}.e\n";
			print JOB "#\$ -o $dta_path/${job_name}_${i}.o\n";			
			foreach (@dta_file_arrays) {
				print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $parameter -dta_path $dta_path $_\n";	
			}
		} else {
			print JOB "perl $dta_path/runsearch_shell.pl -job_num $i -param $parameter -dta_path $dta_path $dta_file_temp\n";	
		}
		close(JOB);
	}


		
	## Run/control parallel jobs 
	my $job_list;
	if ($params -> {'cluster'} eq '1') {	## Cluster system
		if ($params -> {'Job_Management_System'} eq 'LSF') {
			for (my $i = 0; $i < $job_num; $i++) {
				my $command_line = qq(cd $dta_path && bsub <${job_name}_${i}.sh);
				my $job = qx[$command_line];
#				print "return code is ",$?,"\n";
				chomp $job;
				my $job_id = 0;
				if ($job =~ /Job \<(\d*)\> is/) {
					$job_id = $1;
				}
				$job_list -> {$job_id} = 1;
			}
		} elsif ($params -> {'Job_Management_System'} eq 'SGE') {
			for (my $i = 0; $i < $job_num; $i++) {
				my $job_name = "${job_name}_${i}.sh";
				my $command_line = qq(cd $dta_path && qsub -cwd $job_name);
				my $job = qx[$command_line];
				chomp $job;
				my $job_id = 0;
				if ($job =~ /$job_name \<(\d*)\> is/) {
					$job_id = $1;
				}
				$job_list -> {$job_id} = 1;
				my $count = $i + 1;
				print "\r  $count jobs were submitted";
			}
		}
		print "\n  You submitted $job_num jobs for database search\n";
		Check_Job_stat("${job_name}_", $job_num, $dta_path,$job_list);
	} elsif($params->{'cluster'} eq '0') {	## Non-cluster system, but multiple cores
		my $pm = new Parallel::ForkManager($MAX_PROCESSES);
		for my $i (0 .. $MAX_PROCESSES) {
			$pm -> start and next;
			my $job_name = "${job_name}_${i}.sh";			
			system ("cd $dta_path && sh $job_name >/dev/null 2>&1");
			print "\r  $i jobs were submitted";	
			Check_Job_stat("${job_name}_", $job_num, $dta_path);
			$pm -> finish; # Terminates the child process
		}
		$pm -> wait_all_children;		
	}

	## Checking accidentally unfinished jobs  
	my @temp_file_array;
	for (my $k = 0; $k <= $#$file_array; $k++) {
		my $data_file = $file_array -> [$k];
		my $out_file = $data_file;
		$out_file =~ s/\.dta$/\.spout/;
		if (!-e $out_file) {
			push (@temp_file_array, $data_file);
		}
	}
	return \@temp_file_array;
}

sub LaunchParallelJob {
	my ($FileName, $cmd, $GridType, $outputName, $dta_path) = @_;
	open (JOB, ">$FileName") || die "can not open $FileName\n";
	if ($GridType eq 'LSF') {
		print JOB "#BSUB -P prot\n";
		print JOB "#BSUB -q normal\n";
		print JOB "#BSUB -M 20000\n";
		print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
		print JOB "#BSUB -eo $dta_path/$outputName.e\n";
		print JOB "#BSUB -oo $dta_path/$outputName.o\n";
		print JOB $cmd;		
		close (JOB);
		system (qq(bsub <$FileName >/dev/null 2>&1));	
	}
	if ($GridType eq 'SGE') {
		print JOB "#!/bin/bash\n";
		print JOB "#\$ \-S /bin/bash\n";  #In our cluster this line is esential for executing some bash commands such as for
		print JOB "#\$ \-N $outputName\n";
		print JOB "#\$ \-e $dta_path/$outputName.e\n";
		print JOB "#\$ \-o $dta_path/$outputName.o\n";
		print JOB $cmd;
		close (JOB);
		system (qq(qsub -cwd -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G $FileName >/dev/null 2>&1));
	}
	if ($GridType eq 'PBS') {
		print JOB "#!/bin/bash\n";
		print JOB "#PBS -N $outputName\n";
		print JOB "#PBS -e $dta_path/$outputName.e\n"; 
		print JOB "#PBS -o $dta_path/$outputName.o"; 			
		print JOB $cmd;
		close (JOB);
		system (qq(qsub -cwd $FileName >/dev/null 2>&1));
	}
	close(JOB);
}

sub MergeFile {
	my $nbRanges = shift;
	my $cmd = "for i in {0..$nbRanges} do; cat "
}
	
sub Check_Job_stat {
	## Check the status of distributed jobs
	my ($jobs_prefix, $job_num, $dta_path, $jobs_hashref) = @_;
	my $job_info = 1;
    my ($username) = getpwuid($<);
	my $command_line = "";
	my $dot = ".";
	while ($job_info) {
		if ($params -> {'cluster'} eq 0) {
			if ($jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/) {
				my @outfile = glob("$dta_path\/\*.spout");
				my @dtafile = glob("$dta_path\/\*.dta");
				my $outfile_num = scalar(@outfile);
				my $dtafile_num = scalar(@dtafile);
				print "\r  $outfile_num files have done         ";				
				if ($outfile_num == $dtafile_num) {
					$job_info = 0;
				} else {
					$job_info = 1;
				}
				sleep(30);
			}
		}
		if ($params -> {'Job_Management_System'} eq 'LSF') {
			$command_line =  "bjobs -u $username";
		} elsif ($params -> {'Job_Management_System'} eq 'SGE') {
			$command_line =  "qstat -u $username";
		}
		my $job_status = qx[$command_line 2>&1];
		my @job_status_array = split(/\n/, $job_status);
	
		## Consider only the ones that we submitted
		if ($params -> {'Job_Management_System'} eq 'LSF') {
			$command_line =  "bjobs -u $username -noheader";
			my $job_status = qx[$command_line];
			my $running_jobs = set_intersection( $jobs_hashref, parse_bjobs_output($job_status) );
			my $job_number = $job_num - scalar(@$running_jobs);
			if( scalar(@$running_jobs) == 0 ) {
			    print "\r  $job_num jobs finished          ";			    
			    $job_info = 0;
			}
			else {
			    print "\r  $job_number jobs finished          ";			    
			}
			sleep(10);
			
		} elsif ($params -> {'Job_Management_System'} eq 'SGE') {
			@job_status_array = grep(/$jobs_prefix/, @job_status_array);
			if ($job_status =~ /No unfinished job found/) {
				$job_info = 0;
				print "  \n";
			} elsif (scalar(@job_status_array) == 0) {
				$job_info = 0;
			} elsif ($job_status_array[1] =~ /PEND/) {
				print "\r cluster is busy, please be patient!          ";
				sleep(100);
			} elsif ($jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/) {
				my $check_command = "ls -f $dta_path\/\*.spout \| wc -l";
				my @outfile = glob("$dta_path\/\*.spout");
				my $outfile_num = scalar @outfile;
				print "\r  $outfile_num files have done         ";
				sleep(30);
			} elsif ($jobs_prefix eq "job_db_" || $jobs_prefix eq "sort_" || $jobs_prefix eq "merge_"  || $jobs_prefix =~ /sch_/ || $jobs_prefix =~ /resch_/) {
				if ($params -> {'Job_Management_System'} eq 'LSF') {	
					$command_line =  "bjobs -u $username";
				} elsif ($params -> {'Job_Management_System'} eq 'SGE') {
					$command_line =  "qstat -u $username";
				}
				my $job_status = qx[$command_line];
				my @job_status_array = split(/\n/, $job_status);
				my $job_number = $job_num - scalar(@job_status_array) + 2;
				if (scalar(@job_status_array) == 0) {
					print "\r  $job_num jobs finished          ";
				} else {
					print "\r  $job_number jobs finished          ";
				}
				if (scalar(@job_status_array) > 0) {
					$job_info = 1;
				} else {
					print "\r  $job_num jobs finished          ";
					$job_info = 0;
				}			
			}			
		}
	}
}


sub database_creation {
	my $databasename = $params -> {'database_name'};
	my $database_path = dirname($databasename); 
	my $database_basename = basename($databasename);
	my $digestion = $params -> {'digestion'};
	my $num_dynamic_mod = 1;
	my $random = int(rand(100));
	my $tmp_database_path = $database_path . "/.tmp$random/";
	foreach my $key (keys %$params) {
		if ($key =~ /dynamic_/) {
			$num_dynamic_mod++;
		}
	}
	if ($params -> {search_engine} eq 'SEQUEST') {
		## Do nothing?
	} elsif (!(-e($databasename) and $databasename =~ /.mdx/)) {
		print "  Creating database\n";
		my $mass_index = $databasename . ".mdx";
		my $peptide_index = $databasename . ".pdx";
		my $protein_index = $databasename . ".prdx";
		if (-e($mass_index)) {
			print "  Do you want to remove the old database with same name? (y/n): ";
			chomp (my $choice = <STDIN>);
			if ($choice eq "yes" || $choice eq "y") {
				print "  Removing old database\n";
				system (qq(rm -f $mass_index >/dev/null 2>&1));
				system (qq(rm -f $peptide_index >/dev/null 2>&1));
				system (qq(rm -f $protein_index >/dev/null 2>&1));
			}
		}
		$params -> {'database_name'} = $mass_index;

		## Create database files using cluster system
		my $total_protein_number = 0;
		open (FASTA, $databasename) || die "can not open the database\n";
		while (<FASTA>) {
			$total_protein_number++ if($_ =~ /^\>/);
		}
		close (FASTA);
		my $num_protein_per_job = int($total_protein_number / (200 * $num_dynamic_mod)) + 1;
		my $protein_num=0;
		my $k=0;
		open (FASTA, $databasename) || die "can not open the database\n";
		system (qq(mkdir $tmp_database_path >/dev/null 2>&1));
		while (<FASTA>) {
			$protein_num++ if($_ =~ /^\>/);
			if (($protein_num % $num_protein_per_job) == 1) {
				if ($_ =~ /^\>/) {
					$k++;
					print "\r  Generating $k temporary files";
					open (FASTATEMP, ">$tmp_database_path/temp_${k}_${database_basename}");
				}
				print FASTATEMP "$_";
			} else {
				print FASTATEMP "$_";
			}
		}
		close(FASTA);
		print "\n";
		
		## Create job files
		my $job = new Spiders::Job();
		my $abs_parameter = abs_path($parameter);
		$job -> set_library_path($library);			
		my $num_mass_region = 40;
		$job -> make_createdb_script($tmp_database_path, $abs_parameter, $num_mass_region, $num_protein_per_job);

		## Submet the job files
		my $job_num = int($protein_num / $num_protein_per_job) + 1;	
		for (my $i = 1; $i <= $job_num; $i++) {
			open (JOB, ">$tmp_database_path/job_db_$i.sh") || die "can not open the job db files\n";
			if ($params -> {'Job_Management_System'} eq 'LSF') {
				print JOB "#BSUB -P prot\n";
				print JOB "#BSUB -q normal\n";
				print JOB "#BSUB -M 20000\n";
				print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
				print JOB "#BSUB -eo $tmp_database_path/$i.e\n";
				print JOB "#BSUB -oo $tmp_database_path/$i.o\n";
				print JOB "perl $tmp_database_path/create_db.pl $tmp_database_path/temp_${i}_${database_basename}\n";
				print JOB "rm -f $tmp_database_path/temp_${i}_${database_basename}";
				close (JOB);
				system (qq(cd $tmp_database_path && bsub <job_db_$i.sh >/dev/null 2>&1));				
			}
			if ($params -> {'Job_Management_System'} eq 'SGE') {
				print JOB "#!/bin/bash\n";
				print JOB "#\$ \-N job_db_$i\n";
				print JOB "#\$ \-e $tmp_database_path/$i.e\n";
				print JOB "#\$ \-o $tmp_database_path/$i.o\n";
				print JOB "perl $tmp_database_path/create_db.pl $tmp_database_path/temp_${i}_${database_basename}\n";
				print JOB "rm -rf $tmp_database_path/temp_${i}_${database_basename}\n";
				close (JOB);
				system (qq(qsub -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G "$tmp_database_path/job_db_$i.sh" >/dev/null 2>&1));
			}
			close (JOB);
		}
		print "  You submit $job_num jobs for creating index files\n";
		Check_Job_stat("job_db_",$job_num);

		## Merge individual results files into a large file
		my $j = 0;
		my $l = 0;
		my $prev_run = 0;
		print "\n  Sorting indexes\n";
		$job -> make_partialidx_script($tmp_database_path);
		my $prot_job_num = int($num_mass_region / 2);
		my $tmp_mass_index = "$tmp_database_path/$database_basename" . ".mdx";
		my $tmp_peptide_index = "$tmp_database_path/$database_basename" . ".pdx";
		my $tmp_protein_index = "$tmp_database_path/$database_basename" . ".prdx";
		for (my $m = 0; $m < $num_mass_region; $m++) { 		
			my ($mass_hash, $peptide_hash, $protein_hash);	
			## Create bash files
			my $parameters = {'range' => $m,				
							  'GridType' => $params -> {'Job_Management_System'},
							  'JobNum' => $job_num,
							  'database_path' => $tmp_database_path,
							  'dta_path' => $tmp_database_path,
							  'mass_index' => $tmp_mass_index . ".$m",
							  'peptide_index' => $tmp_peptide_index . ".$m",
							  'protein_index' => $protein_index,
							   'databasename' => $database_basename,
							   'num_pro_per_job' => $num_protein_per_job,
							   'prot_job_num' => $prot_job_num,
							   'digestion' => $digestion,
							   };
			Create_Sort_BashFile($parameters, $tmp_database_path);							
		}
		print "  You submit $num_mass_region jobs for sorting index files\n";
		Check_Job_stat("sort_",$num_mass_region);
		print "  $num_mass_region jobs finished";		
		print "\n  Mergering files\n";
		my $merge_job_num = $num_mass_region - 1;
		
		my $cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_mass_index.".'$i'." >> $mass_index\n done\n";
		my $FileName = "$tmp_database_path/merge_mass_index.sh";
		LaunchParallelJob($FileName, $cmd, $params -> {'Job_Management_System'}, "merge_mass_index", $tmp_database_path);		
		$cmd = "for i in {0..$merge_job_num} \n do\n cat $tmp_peptide_index.".'$i'." >> $peptide_index\n done\n";
		$FileName = "$tmp_database_path/merge_peptide_index.sh";
		LaunchParallelJob($FileName, $cmd, $params -> {'Job_Management_System'}, "merge_peptide_index", $tmp_database_path);
		
		## Create a database summary file	
		my @summaryfiles = ();
		for (my $i = 0; $i < $merge_job_num; $i++) {
			my $sumfile = "$tmp_mass_index.${i}.summary";
			push (@summaryfiles, $sumfile);
		}
		my $outfile = $mass_index;
		$outfile =~ s/mdx/sdx/;
		summary_db(\@summaryfiles, $outfile);

		print "  You submit 2 jobs for merging index files\n";
		Check_Job_stat("merge_", "2");
		print "  2 jobs finished";
		print "\n  Removing temporary files\n";
		system(qq(rm -rf $tmp_database_path >/dev/null 2>&1));
		print "  Database creation completed\n";
	}
}

sub summary_db {
	my ($files, $outfile) = @_;
	my $fasta_file = $params -> {'database_name'};
	my $enzyme_info = $params -> {'enzyme_info'};
	my $digestion = $params -> {'digestion'};	
	my $mis_cleavage = $params -> {'max_mis_cleavage'};
	my $min_mass = $params -> {'min_peptide_mass'};	
	my $max_mass = $params -> {'max_peptide_mass'};	
	my %hash;
	foreach my $file (@$files) {
		open (FILE, $file);
		while (<FILE>) {
			my ($key, $value) = split(/\:/, $_);
			if (!defined($hash{$key})) {
				$hash{$key} = $value;
			} else {
				$hash{$key} += $value;
			}
		}
	}
	close(FILE);
	open(OUTPUT,">$outfile");
	print OUTPUT "JUMP Database index file\n";
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon, $mday, $hour, $min, $sec;	
	print OUTPUT "Database_name = ", $fasta_file, "\n";
	print OUTPUT "enzyme_info = ", $enzyme_info, "\n";
	print OUTPUT "digestion = ", $digestion, "\n";
	print OUTPUT "max_mis_cleavage = ", $mis_cleavage, "\n";
	print OUTPUT "min_peptide_mass = ", $min_mass, "\n";
	print OUTPUT "max_peptide_mass = ", $max_mass, "\n";
	foreach (sort keys %hash) {
		print OUTPUT $_," = ", $hash{$_}, "\n";
	}
	close (OUTPUT);
}

sub topMS2scans {
	my ($msmshash) = @_;
	my (%hash, @array, $preNumber);
	foreach my $number (keys %{$msmshash}) {
		$preNumber = $$msmshash{$number}{'survey'};
		if (defined($preNumber)) {
			if (defined($hash{$preNumber})) {
				$hash{$preNumber}++;
			} else {
				$hash{$preNumber} = 1;
			}
		}
	}
	for (my $i = 0; $i <= 10; $i++) {
		$array[$i] = 0;
	}
	foreach my $number (keys %hash) {
		for (my $i = 0; $i < $hash{$number}; $i++) {
			$array[$i]++;
		}
	}
	return (\@array);
}

sub delete_run_dir {
	my ($run) = @_;
	my $tmp = int(rand(100000));
	$tmp = "." . $tmp;
	system (qq(mkdir $tmp;cp $run\/*.params $tmp;cp $run\/*.pl $tmp));
	system (qq(rm -rf $run > /dev/null 2>&1));
	system (qq(mkdir $run > /dev/null 2>&1));
	system (qq(mv $tmp/* $run/));
	system (qq(rm -rf $tmp > /dev/null 2>&1));
}

sub usage {
    my $msg = shift;
print <<EOF;
################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump search using JUMP                  ****     # 
#       ****  Version 12.1.0                          ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################

Usage: $progname -p parameterfile rawfile.raw 
	or
       $progname -p parameterfile rawfile.mzXML
	

EOF
if( defined($msg) ) { print STDERR $msg; }
exit 1;
}

sub parse_bjobs_output {
    my $output = shift;
    my @arr = ($output =~ /^(\d+).*$/mg);
    my %rv;
    @rv{@arr} = @arr;
    return \%rv;
}

sub set_intersection {
    my ($s1,$s2) = @_;
    my @rv;
    foreach my $i1 (keys(%$s1)) {
	if(defined($s2->{$i1})) {
	    push( @rv, $i1 );
	}
    }
    return \@rv;
}

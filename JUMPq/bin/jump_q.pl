#!/bin/env perl 

our $VERSION = 1.13.001;
use lib $ENV{"JUMP_Q_LIB"};
#
use File::Basename;
use Cwd;
use File::Spec;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";
use Utils::Parse;
print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump quantifiction                      ****     #
#       ****  Version 1.13.002                        ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF
    unless( scalar(@ARGV) > 0 ) { help(); }

my $queue;
my $mem;
my $dispatch;
GetOptions('--queue=s'=>\$queue, '--mem=s'=>\$mem, '--dispatch=s'=>\$dispatch);

unless(defined($dispatch)) {
    $dispatch = "batch-interactive";
}

#$mem = 8192;
unless(defined($mem)) {
	my $require_file = "/hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.003/JUMP/bin/batch/BatchLib.pl";
	require $require_file;
    my $jump_qj_fullfile = $ARGV[0];
    my $txtsize = get_txtsize_from_jump_qj($jump_qj_fullfile);
    my $nMEM1 = calculate_mem($txtsize,2); # choice: 2,q
    $mem = $nMEM1*1024;
	unless(defined($queue)) {
		$queue = "standard";
		if ($nMEM1>=200) {
			$queue = "large_mem";
		}
	}
}

unless(defined($queue)) {
    $queue = "standard";
}

if( $dispatch eq "batch-interactive" ) {
	my $nMEM1 = int($mem/1024+0.5);
	my $hint1 = "Applying ".$nMEM1." GB RAM in queue <".$queue."> (please be patient)\n";
	print $hint1;
}

my $cmd;
unless(defined($dispatch) && $dispatch eq 'localhost') {
    $cmd="bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_q.pl" . " " . $ARGV[0];
}
else {
    $cmd="_jump_q.pl " . $ARGV[0];
}
system($cmd);

sub help {
	my ($value) = @_;
	if ($value == 0){
		print "\n";
		print "     Usage: jump_q.pl jump_q.params \n";
		print "\n";
	}
	exit;
}
my ($paramFile) = @ARGV;
my $parse = Utils::Parse -> new();
my %params;
$parse -> parseParams($paramFile, \%params);
$params{'save_dir'} = 'quan_' . $params{'save_dir'};
my $saveDir = getcwd() . "\/$params{'save_dir'}";
my $publicationDir = $saveDir."/publications";
my $input_file = $publicationDir."/id_uni_prot_quan.txt";
if (-e $input_file) {
    my @cmd = ("python $Bin/protmat2genemat.py $input_file");
    system(@cmd);}
else {
    print "NOTE: The protein to gene matrix conversion program cannot operate. This could be because the quantification was done at the peptide level. Thus, protein quantification file is absent.\n\n";}    
    


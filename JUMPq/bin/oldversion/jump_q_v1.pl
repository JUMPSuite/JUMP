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

if(!defined($queue) && !defined($mem)) {
    $queue = 'standard';
    $mem = 8192;
}
elsif(!defined($queue) && defined($mem)) { 
    print "\t--mem cannot be used without --queue\n";
    exit(1);
}
elsif(!defined($mem)) {
    $mem = 8192;
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
    


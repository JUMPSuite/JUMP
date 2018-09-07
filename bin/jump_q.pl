#!/bin/env perl 

our $VERSION = 1.13.0;

use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump quantifiction                      ****     #
#       ****  Version 1.13.0                          ****     #
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
GetOptions('--queue=s'=>\$queue, '--memory=s'=>\$mem, '--dispatch'=>\$dispatch);

if(!defined($queue) && !defined($mem)) {
    $queue = 'normal';
    $mem = 200000;
}
elsif(!defined($queue) && defined($mem)) { 
    print "\t--mem cannot be used without --queue\n";
    exit(1);
}
elsif(!defined($mem)) {
    $mem = 200000;
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

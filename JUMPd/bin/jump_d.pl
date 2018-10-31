#!/bin/env perl

use strict;
use Getopt::Long;

my $dispatch;
my $queue;
GetOptions('--dispatch=s'=>\$dispatch,
	   '--queue=s'=>\$queue);

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}

unless(defined($queue)) {
    $queue = "normal";
}

my $cmd;
unless(defined($dispatch) && $dispatch eq 'localhost') {
    $cmd = "bsub -P prot -q $queue -R \"rusage[mem=20000]\" -Ip _jump_d.pl " . join(" ",@ARGV);
}
else {
    $cmd="_jump_d.pl " . $ARGV[0];
}
system($cmd); 

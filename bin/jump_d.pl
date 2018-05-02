#!/bin/env perl

use strict;

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}

my $cmd = "bsub -P prot -q normal -R \"rusage[mem=20000]\" -Ip _jump_d.pl " . join(" ",@ARGV);
system($cmd); 

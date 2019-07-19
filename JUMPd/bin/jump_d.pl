#!/bin/env perl

use strict;
use Getopt::Long;
use Spiders::Config;
use File::Spec;

my $dispatch;
my $config = new Spiders::Config();
GetOptions('--dispatch=s'=>\$dispatch);

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}

my $cmd;
unless(defined($dispatch) && $dispatch eq 'localhost') {
    my $cmd = "bsub -P prot -q " . $config->get("normal_queue") . " -R \"rusage[mem=20000]\" -Ip _jump_d.pl " . join(" ",@ARGV);
}
else {
    $cmd="_jump_d.pl " . $ARGV[0];
}
system($cmd); 

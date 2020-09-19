#!/bin/env perl

sub usage {
print <<"EOF";
Usage: $progname -p parameterfile
	or
       $progname -p parameterfile
EOF
exit 1;
}

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump localization                       ****     #
#       ****  Version 1.13.002                        ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF

use Getopt::Long;

my ($help,$parameter,$raw_file);
my $queue;
my $mem;
GetOptions('-help|h'=>\$help,
	   '-p=s'=>\$parameter,
	   '--queue=s'=>\$queue, '--mem=s'=>\$mem);

if(!defined($queue) && !defined($mem)) {
    $queue = 'standard';
    $mem = 200000;
}
elsif(!defined($queue) && defined($mem)) { 
    print "\t--mem cannot be used without --queue\n";
    exit(1);
}
elsif(!defined($mem)) {
    $mem = 200000;
}

usage() if ($help || !defined($parameter));

my $cmd="bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_l.pl" . " -p " . $parameter;
system($cmd);

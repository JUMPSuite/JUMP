#!/usr/bin/perl 

our $VERSION = 12.1.0;

use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump filter                             ****     #
#       ****  Version 1.13.0                          ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF
my $queue;
my $mem;
GetOptions('--queue=s'=>\$queue, '--memory=s'=>\$mem);

if(defined($queue) && !defined($mem)) {
    $mem = 262144;
}
elsif(!defined($queue)) { 
    $queue = 'large_mem';
    $mem = 2097152;
}

unless( scalar(@ARGV) > 0 ) { print "\tusage: jump_f.pl <parameter file>\n"; exit(1); }
$cmd="bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_f.pl" . " " . $ARGV[0];
system($cmd);

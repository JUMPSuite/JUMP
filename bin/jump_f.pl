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

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
my $outname = sprintf("%4d-%02d-%02d_%02d:%02d:%02d_",$year+1900,$mon+1,$mday,$hour,$min,$sec);	
$outname .= $ARGV[0] . ".out";

unless( scalar(@ARGV) > 0 ) { print "\tusage: jump_f.pl <parameter file>\n"; exit(1); }
$cmd="bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_f.pl" . " " . $ARGV[0] . " 2>&1 | tee $outname";
system($cmd);

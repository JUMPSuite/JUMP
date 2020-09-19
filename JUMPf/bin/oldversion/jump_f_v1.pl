#!/hpcf/apps/perl/install/5.10.1/bin/perl 

eval 'exec /hpcf/apps/perl/install/5.10.1/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

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
#       ****  Version 1.13.001                        ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF
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

unless( scalar(@ARGV) > 0 ) { print "\tusage: jump_f.pl <parameter file>\n"; exit(1); }

my $cmd;
unless(defined($dispatch) && $dispatch eq 'localhost') {
    $cmd="bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_f.pl" . " " . $ARGV[0];
}
else {
    $cmd="_jump_f.pl " . $ARGV[0];
}
system($cmd);

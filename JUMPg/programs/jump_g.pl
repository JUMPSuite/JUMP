our $VERSION = 2.3.5;

use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;
use Spiders::Config;
my $config = new Spiders::Config();

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump G                                  ****     #
#       ****  Version 1.13.1                          ****     #
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
    $queue = $config->get("normal_queue");
    $mem = 200000;
}
elsif(!defined($queue) && defined($mem)) { 
    print "\t--mem cannot be used without --queue\n";
    exit(1);
}
elsif(!defined($mem)) {
    $mem = 200000;
}

unless( scalar(@ARGV) > 0 ) { print "\tUsage: perl rundtas rundtas.params qc_MSMS_input.txt\n"; }
$cmd="bsub -env all d -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip JUMPg_v2.3.pl " . join( " ", @ARGV );
system($cmd);

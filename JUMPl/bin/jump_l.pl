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
#       ****  Version 1.13.1                          ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF

use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;

my ($help,$parameter,$raw_file);
my $dispatch;
GetOptions('-help|h'=>\$help,
	   '--dispatch'=>\$dispatch );

usage() if ($help || !defined($parameter));

my $cmd;
if((defined($dispatch) && $dispatch eq "localhost") || Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_LOCALIZATION);
    $cmd="$batchCmd _jump_l.pl" . " " . $ARGV[0];
} elsif(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="_jump_l.pl " . $ARGV[0];
}
system($cmd);

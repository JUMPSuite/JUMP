#!/usr/bin/perl 

our $VERSION = 12.1.0;

use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump filter                             ****     #
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
my $dispatch;
my $config = new Spiders::Config();
GetOptions('--queue=s'=>\$queue, '--memory=s'=>\$mem, '--dispatch=s'=>\$dispatch);

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

unless( scalar(@ARGV) > 0 ) { print "\tusage: jump_f.pl <parameter file>\n"; exit(1); }

my $cmd;
if(defined($dispatch) || Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_FILTER);
    $cmd="$batchCmd _jump_f.pl" . " " . $ARGV[0];
} elsif(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="_jump_f.pl " . $ARGV[0];
}
system($cmd);

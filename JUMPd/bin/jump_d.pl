#!/bin/env perl

use Carp;
use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
use File::Spec;
use Spiders::Which;

my $dispatch;
my $config = new Spiders::Config();
GetOptions('--dispatch=s'=>\$dispatch);

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}

my $cmd;
my $jumpd = Spiders::Which::which( "_jump_d.pl" );
if (Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_DATABASE);    
    $cmd = "$batchCmd \"perl $jumpd" . " " . $ARGV[0] ."\"";
} elsif ((defined($dispatch) && $dispatch eq "localhost") || Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->SMP) {
    $cmd = "perl $jumpd " . $ARGV[0];
}
0 == system($cmd) || croak("command \"$cmd\" failed to execute with code $?");

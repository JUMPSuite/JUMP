#!/bin/env perl

use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
use File::Spec;

my $dispatch;
my $config = new Spiders::Config();
GetOptions('--dispatch=s'=>\$dispatch);

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}

my $cmd;
if ((defined($dispatch) && $dispatch eq "localhost") || Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_DATABASE);    
    $cmd = "$batchCmd \"_jump_d.pl" . " " . $ARGV[0] ."\"";
} elsif (Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->SMP) {
    $cmd = "_jump_d.pl " . $ARGV[0];
}
system($cmd); 

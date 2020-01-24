#!/bin/env perl

use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
use File::Spec;

my $queue;
my $mem;
my $dispatch;
my $config = new Spiders::Config();
GetOptions('--dispatch=s'=>\$dispatch, '--queue=s'=>\$queue, '--memory=s'=>\$mem);

if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}
if (!defined($queue) && !defined($mem)) {
	$queue = $config->get("normal_queue");
	$mem = 200000;
} elsif (!defined($queue) && defined($mem)) {
	print "\t--mem cannot be used without --queue\n";
	exit(1);
} elsif (!defined($mem)) {
	$mem = 200000;
}

my $cmd;
if (defined($dispatch) || Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_DATABASE);    
    $cmd = "$batchCmd \"_jump_d.pl" . " " . $ARGV[0] ."\"";
} elsif (Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->SMP) {
    $cmd = "_jump_d.pl " . $ARGV[0];
}
system($cmd); 

#!/bin/env perl

use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
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
	$cmd = "bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_d.pl" . " " . $ARGV[0];
} elsif (Spiders::ClusterConfig::getClusterConfig($config, $params) eq Spiders::ClusterConfig->SMP) {
	$cmd = "_jump_d.pl " . $ARGV[0];
}
system($cmd); 

our $VERSION = 2.3.5;

use Carp;
use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;
use Spiders::Config;
use Spiders::Which;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;

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
my $dispatch;
my $config = new Spiders::Config();
GetOptions('--dispatch=s'=>\$dispatch);

unless( scalar(@ARGV) > 0 ) { print "\tUsage: perl rundtas rundtas.params qc_MSMS_input.txt\n"; }
my $cmd="JUMPg_v2.3.pl " . join( " ", @ARGV );
my $jumpg = Spiders::Which::which( "JUMPg_v2.3.pl" );
if(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_G);
    $cmd="$batchCmd perl $jumpg ";
} elsif((defined($dispatch) && $dispatch eq "localhost") ||
	Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="perl $jumpg ";
}
0 == system($cmd) || croak("command \"$cmd\" failed to execute with code $?");

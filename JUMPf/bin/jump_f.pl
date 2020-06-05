our $VERSION = 12.1.0;

use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
use Spiders::Which;
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
my $dispatch;
my $config = new Spiders::Config();
GetOptions('--dispatch=s'=>\$dispatch);

unless( scalar(@ARGV) > 0 ) { print "\tusage: jump_f.pl <parameter file>\n"; exit(1); }

my $cmd;
my $jumpf = Spiders::Which::which( "_jump_f.pl" );
if(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_FILTER);
    $cmd="$batchCmd perl $jumpf " . $ARGV[0];
} elsif((defined($dispatch) && $dispatch eq "localhost") ||
	Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="perl $jumpf " . $ARGV[0];
}
system($cmd);

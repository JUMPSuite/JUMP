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
#       ****  Copyright (C) 2012 - 2020               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF

use Carp;
use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
use Spiders::Which;
my $config = new Spiders::Config();

my ($help,$parameter,$raw_file);
my $dispatch;
GetOptions('-help|h'=>\$help,
	   '-p=s'=>\$parameter,
	   '--dispatch=s'=>\$dispatch );

usage() if ($help || !defined($parameter));

my $cmd;
my $jumpl = Spiders::Which::which("_jump_l.pl");
if((defined($dispatch) && $dispatch eq "localhost") || Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_LOCALIZATION);
    $cmd="$batchCmd perl $jumpl" . " -p $parameter"
} elsif(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="perl $jumpl " . " -p $parameter"
}
system($cmd) || croak("command \"$cmd\" failed to execute with code $?");

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
use Spiders::Which;

my ($help,$parameter,$raw_file);
my $dispatch;
GetOptions('-help|h'=>\$help,
	   '--dispatch'=>\$dispatch );

usage() if ($help || !defined($parameter));

my $cmd;
my $jumpl = Spiders::Which::which("jump_l.pl");
if((defined($dispatch) && $dispatch eq "localhost") || Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_LOCALIZATION);
    $cmd="$batchCmd $jumpl" . " " . $ARGV[0];
} elsif(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="$jumpl " . $ARGV[0];
}
system($cmd);

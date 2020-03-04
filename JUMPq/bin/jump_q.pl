our $VERSION = 1.13.1;

use lib $ENV{"JUMP_Q_LIB"};
use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::BatchSystem;
use Spiders::Which;
use Utils::Parse;
my $config = new Spiders::Config();

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump quantifiction                      ****     #
#       ****  Version 1.13.1                          ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF
    unless( scalar(@ARGV) > 0 ) { help(); }

my $dispatch;
GetOptions('--dispatch=s'=>\$dispatch);

my %params;
Utils::Parse->new($ARGV[0],\%params);

my $cmd;
my $jumpq = Sipders::Which::which("_jump_q.pl");
if((defined($dispatch) && $dispatch eq "localhost") || Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_QUANTIFICATION);
    $cmd="$batchCmd perl $jumpq " . $ARGV[0];
} elsif(Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $cmd="perl $jumpq " . $ARGV[0];
}
system($cmd);

sub help {
	my ($value) = @_;
	if ($value == 0){
		print "\n";
		print "     Usage: jump_q.pl jump_q.params \n";
		print "\n";
	}
	exit;
}

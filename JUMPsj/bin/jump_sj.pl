#!/bin/env perl

my $Bin=$ENV{"JUMP_SJ_LIB"};
use lib $ENV{"JUMP_SJ_LIB"};
use Getopt::Long;
use Spiders::JUMPmain;
use File::Temp;
use Cwd;
use Cwd 'abs_path';
use Spiders::Config;
use Spiders::ClusterConfig;
use Spiders::Params;
use Spiders::JobManager;
my $config = new Spiders::Config();
our $VERSION = 1.13.1;

my ($help,$parameter,$raw_file,$dispatch);
my %options;
GetOptions('-help|h'=>\$help, '--dispatch=s'=>\$dispatch,
	   '-p=s'=>\$parameter, '--dtafile-location=s'=>\${$options{'--dtafile-location'}},
	   '--keep-dtafiles'=>\${$options{'--keep-dtafiles'}},
	   '--queue=s'=>\${$options{'--queue'}},
	   '--preserve-input'=>\${$options{'--preserve-input'}},
	   '--max-jobs=s'=>\${$options{'--max-jobs'}}
    );

my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();

unless(defined(${$options{'--dtas-backend'}})) {
    ${$options{'--dtas-backend'}} = 'idx';
}

if(!defined($dispatch) && Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->CLUSTER) {
    $dispatch = ($config->get('compute_on_login_node') eq '1' ? "localhost" : 'batch-interactive');
} elsif(!defined($dispatch) && Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP) {
    $dispatch = "localhost";
}

unless(defined(${$options{'--max-jobs'}})) {
    ${$options{'--max-jobs'}} = 512;
}


unless(defined(${$options{'--queue'}})) {
    ${$options{'--queue'}} = $config->get("normal_queue");
}
my $queue = ${$options{'--queue'}};


if(defined(${$options{'--dtafile-location'}}) && !File::Spec->file_name_is_absolute(${$options{'--dtafile-location'}})) {
    ${$options{'--dtafile-location'}} = File::Spec->rel2abs(${$options{'--dtafile-location'}});
}

usage() if ($help || !defined($parameter));
#check_input($raw_file,\$parameter);

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump search using JUMP                  ****     #
#       ****  Version 1.13.1                          ****     #
#       ****  Xusheng Wang / Junmin Peng              ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################


EOF

for my $k (keys(%options)) {
    if(defined(${$options{$k}})) {
	$options_str .= $k . "=" . ${$options{$k}} . " ";
    }
}

if( $dispatch eq "batch-interactive" && $config->get('compute_on_login_node') eq '0' ) {
    my ($handle,$tname) = File::Temp::mkstemp( "JUMPSJXXXXXXXXXXXXXX" );
    my $batchSystem = new Spiders::BatchSystem();
    my $batchCmd = $batchSystem->getBatchCmd(Spiders::BatchSystem->JUMP_SEARCH);
    my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter . " " . $options_str;
    system( "$batchCmd \"$cmd --dispatch=localhost 2>&1 | tee $tname ; jump_sj_log.pl < $tname ; rm $tname\"" );
}
# elsif( $dispatch eq "batch" ) {
#     my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter . " " . $options_str;
#     system( "bsub -env all -P prot -q $queue -R \"rusage[mem=32768]\" \"$cmd --dispatch=localhost | jump_sj_log.pl\"" );
# }
elsif( $dispatch eq "localhost" || 
       Spiders::ClusterConfig::getClusterConfig($config,$params) eq Spiders::ClusterConfig->SMP || 
       $config->get('compute_on_login_node') eq '1' ) {
    my $library = $Bin;
    my $main = new Spiders::JUMPmain();
    $main->set_library($library);
    $main->main($parameter,\@ARGV,\%options);
}
else {
    print "argument to --dispatch must be one of \"batch-interactive\", or \"localhost\"\n";
    exit -1;
}
sub usage {

print <<"EOF";
	
################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump search using JUMP                  ****     #
#       ****  Version 1.13.1                          ****     #
#       ****  Xusheng Wang / Junmin Peng              ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################


Usage: $progname -p parameterfile rawfile.raw 
	or
       $progname -p parameterfile rawfile.mzXML
	

EOF
exit 1;
}

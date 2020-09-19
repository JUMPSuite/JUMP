#!/bin/env perl

my $Bin=$ENV{"JUMP_SJ_LIB"};
use lib $ENV{"JUMP_SJ_LIB"};
use Getopt::Long;
use Spiders::JUMPmain2;
use File::Temp;
use Cwd;
use Cwd 'abs_path';
our $VERSION = 1.13.0;

my ($help,$parameter,$raw_file,$dispatch,$mem);
my %options;
GetOptions('-help|h'=>\$help, '--dispatch=s'=>\$dispatch,
	   '-p=s'=>\$parameter, '--dtafile-location=s'=>\${$options{'--dtafile-location'}},
	   '--keep-dtafiles'=>\${$options{'--keep-dtafiles'}},
	   '--queue=s'=>\${$options{'--queue'}},
       '--mem=s'=>\${$options{'--mem'}},
	   '--preserve-input'=>\${$options{'--preserve-input'}},
	   '--max-jobs=s'=>\${$options{'--max-jobs'}}
    );

unless(defined($dispatch)) {
    #$dispatch = "batch-interactive";
    $dispatch = "localhost";
}

unless(defined(${$options{'--max-jobs'}})) {
    ${$options{'--max-jobs'}} = 512;
}


unless(defined(${$options{'--queue'}})) {
    ${$options{'--queue'}} = "standard";
}
my $queue = ${$options{'--queue'}};

#$mem = 4096;
unless(defined(${$options{'--mem'}})) {
    ${$options{'--mem'}} = 12288;
}
$mem = ${$options{'--mem'}};


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
#       ****  Version 1.13.0                          ****     #
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

if( $dispatch eq "batch-interactive" ) {
    my ($handle,$tname) = File::Temp::mkstemp( "JUMPSJXXXXXXXXXXXXXX" );
    my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter . " " . $options_str;
    system( "bsub -env all -P prot -q $queue -R \"rusage[mem=$mem]\" -Is \"$cmd --dispatch=localhost 2>&1 | tee $tname ; jump_sj_log.pl < $tname ; rm $tname\"" );
}
elsif( $dispatch eq "batch" ) {
    my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter . " " . $options_str;
    system( "bsub -env all -P prot -q $queue -R \"rusage[mem=$mem]\" \"$cmd --dispatch=localhost | jump_sj_log.pl\"" );
}
elsif( $dispatch eq "localhost" ) {
    my $library = $Bin;
    my $main = new Spiders::JUMPmain2();
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
#       ****  Version 1.13.0                          ****     #
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

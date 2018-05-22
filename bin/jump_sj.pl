#!/bin/env perl

my $Bin=$ENV{"JUMP_SJ_LIB"};
use lib $ENV{"JUMP_SJ_LIB"};
use Getopt::Long;
use strict;
use Spiders::JUMPmain;
use File::Temp;
use Cwd;
use Cwd 'abs_path';
our $VERSION = 1.13.0;

my ($help,$parameter,$raw_file,$dispatch,$queue);
my %options;
GetOptions('-help|h'=>\$help, 
	   '--dispatch=s'=>\$dispatch,
	   '-p=s'=>\$parameter, 
	   '--dtafile-location=s'=>\${$options{'--dtafile-location'}},
	   '--keep-dtafiles'=>\${$options{'--keep-dtafiles'}}, 
	   '--queue=s'=>\$queue, 
	   '--dtas-backend=s'=>\${$options{'--dtas-backend'}},
	   '--preserve-input'=>\${$options{'--preserve-input'}}
    );

unless(defined(${$options{'--dtas-backend'}})) {
    ${$options{'--dtas-backend'}} = 'idx';
}

unless(defined($queue)) {
    $queue = 'heavy_io';
}

unless(defined($dispatch)) {
    $dispatch = "batch-interactive";
}

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

my $options_str;
for my $k (keys(%options)) {
    if(defined(${$options{$k}})) {
	$options_str .= $k . "=" . ${$options{$k}} . " ";
    }
}

if( $dispatch eq "batch-interactive" ) {
    my ($handle,$tname) = File::Temp::mkstemp( "JUMPSJXXXXXXXXXXXXXX" );
    my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter . " " . $options_str;
    system( "bsub -env all -P prot -q normal -R \"rusage[mem=32768]\" -Is \"$cmd --dispatch=localhost 2>&1 | tee $tname ; jump_sj_log.pl < $tname ; rm $tname\"" );
}
elsif( $dispatch eq "batch-parallel" ) {
    # foreach my $arg (@ARGV) {
    # 	my $cmd = 'jump_sj.pl ' . $arg . " -p " . $parameter . " " . $options_str;
    # 	print "submitting job for $arg\n";
    # 	system( "$cmd --dispatch=batch --queue=$queue &> /dev/null" );
    # }
    print "batch-parallel mode not yet supported\n";
    exit -1;
}
elsif( $dispatch eq "batch" ) {
    my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter . " " . $options_str;
    system( "bsub -env all -P prot -q normal -R \"rusage[mem=32768]\" \"$cmd --dispatch=localhost | jump_sj_log.pl\"" );
}
elsif( $dispatch eq "localhost" ) {
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
#       ****  Version 1.13.0                          ****     #
#       ****  Xusheng Wang / Junmin Peng              ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################


Usage: $0 -p parameterfile rawfile.raw 
	or
       $0 -p parameterfile rawfile.mzXML
	

EOF
exit 1;
}

#!/bin/env perl

my $Bin=$ENV{"JUMP_SJ_LIB"};
use lib $ENV{"JUMP_SJ_LIB"};
use Getopt::Long;
use Spiders::JUMPmain;
use Cwd;
use Cwd 'abs_path';
our $VERSION = 1.13.0;

my ($help,$parameter,$raw_file,$nobatch);
GetOptions('-help|h'=>\$help, '--no-batch'=>\$nobatch,
			'-p=s'=>\$parameter,
		);

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

unless( defined($nobatch) ) {
    my $cmd = 'jump_sj.pl ' . join( ' ', @ARGV ) . " -p " . $parameter;
    print "bsub -env all -P prot -q normal -R \"rusage[mem=32768]\" -Is $cmd --no-batch",'\n';
    system( "bsub -env all -P prot -q normal -R \"rusage[mem=32768]\" -Is $cmd --no-batch" );
}
else {
    my @args;
    foreach my $arg (@ARGV) {
	unless( $arg =~ /--no-batch/ ) {
	    push( @args, $arg );
	}
    }
    my $library = $Bin;
    my $main = new Spiders::JUMPmain();
    $main->set_library($library);
    $main->main($parameter,\@args);
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

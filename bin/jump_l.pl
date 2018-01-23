#!/bin/env perl

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
#       ****  jump quantifiction                      ****     #
#       ****  Version 1.13.0                          ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF

use Getopt::Long;

my ($help,$parameter,$raw_file);
GetOptions('-help|h'=>\$help,
	'-p=s'=>\$parameter,
);

usage() if ($help || !defined($parameter));

my $cmd="bsub -P prot -q large_mem -R \"rusage[mem=2097152]\" -Ip _jump_l.pl" . " -p " . $parameter;
system($cmd);

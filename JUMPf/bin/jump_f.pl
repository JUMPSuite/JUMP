#!/hpcf/apps/perl/install/5.10.1/bin/perl 

eval 'exec /hpcf/apps/perl/install/5.10.1/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

our $VERSION = 12.1.0;

use File::Basename;
use Cwd 'abs_path';
use File::Spec;
use Getopt::Long;

print <<EOF;

################################################################
#                                                              #
#       **************************************************     #
#       ****                                          ****     #
#       ****  jump filter                             ****     #
#       ****  Version 1.13.001                        ****     #
#       ****  Copyright (C) 2012 - 2017               ****     #
#       ****  All rights reserved                     ****     #
#       ****                                          ****     #
#       **************************************************     #
#                                                              #
################################################################
EOF
my $queue;
my $mem;
my $dispatch;
GetOptions('--queue=s'=>\$queue, '--mem=s'=>\$mem, '--dispatch=s'=>\$dispatch);

if(!defined($queue) && !defined($mem)) {
    $queue = 'standard';
    $mem = 8192;
}
elsif(!defined($queue) && defined($mem)) { 
    print "\t--mem cannot be used without --queue\n";
    exit(1);
}
elsif(!defined($mem)) {
    $mem = 8192;
}

unless( scalar(@ARGV) > 0 ) { print "\tusage: jump_f.pl <parameter file>\n"; exit(1); }

if (!defined($dispatch)) {
    $dispatch = "automem";
}
my $require_file = "/hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.003/JUMP/bin/batch/BatchLib.pl";
require $require_file;
if (defined($dispatch) && $dispatch eq 'automem') {
    my $jump_fj_fullfile = $ARGV[0];
    my $engine = get_engine_from_jump_fj($jump_fj_fullfile);
    my $pepXMLsize = get_pepXMLsize_from_jump_fj($jump_fj_fullfile,$engine);
    my $nMEM1 = calculate_mem($pepXMLsize,1); # choice: 1,f
    $queue = "standard";
    if ($nMEM1>=200) {
        $queue = "rhel7_large_mem";
    }
    $mem = $nMEM1*1000;
    my $hint1 = "Applying ".$nMEM1." GB RAM in queue <".$queue."> (please be patient)\n";
    print $hint1;
}

my $cmd;
unless(defined($dispatch) && $dispatch eq 'localhost') {
    $cmd="bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip _jump_f.pl" . " " . $ARGV[0];
}
else {
    $cmd="_jump_f.pl " . $ARGV[0];
}
system($cmd);

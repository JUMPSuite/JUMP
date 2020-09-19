#!/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin"; # $Bin = '/hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.003/JUMP/bin'


my $require_file = $Bin."/batch/BatchLib.pl"; # use the file $Bin/batch/BatchLib.pl
require $require_file;

# check input
if (scalar(@ARGV)!=1)
{
	die "Usage: jump -batch-id jump_batchID.params\n";
}

# current datapath (e.g. /home/zyuan1/testBatch)
my $cwd_dir = getcwd();

# get nbatch from batchID, and check if they are the same
my $batchid_param_fullfile = $cwd_dir."/".$ARGV[0];
my ($nbatch_ID,$nbatch_fj_params) = get_nbatch_batchID($batchid_param_fullfile);
if ($nbatch_ID!=$nbatch_fj_params) {
	die "In ".$ARGV[0].", some batches in 'input_path_batch' do not have filter results.\n";
}

# write_step1_fbatch
my $fbatch_sh_fullfile = $cwd_dir."/ParameterFiles/Batch/run1_fbatch.sh";
write_step1_fbatch($fbatch_sh_fullfile,$cwd_dir,$Bin,$ARGV[0]);

# first get mem and queue, and then run step1_fbatch
my $pepXMLsize = get_pepXMLsize_from_fbatch($batchid_param_fullfile);
my $nMEM1 = calculate_mem($pepXMLsize,4); # choice: 4,fbatch
submit_jobs($nMEM1,$fbatch_sh_fullfile);
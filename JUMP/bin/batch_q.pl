#!/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin"; # $Bin = '/hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.003/JUMP/bin'
use File::Path qw(make_path remove_tree);

my $require_file = $Bin."/batch/BatchLib.pl"; # use the file $Bin/batch/BatchLib.pl
require $require_file;

# check input
if (scalar(@ARGV)!=1)
{
	die "Usage: jump -batch-q jump_batchQ.params\n";
}

# current datapath (e.g. /home/zyuan1/testBatch)
my $cwd_dir = getcwd();

# get nbatch from batchID and batchQ, and check if they are the same
my $batchq_param_fullfile = $cwd_dir."/".$ARGV[0];
my ($nbatch_Q,$nbatch_internal,$input_mode,$batchID_fullfile,$output_folder_Q,$first1,$second1) = get_nbatch_batchQ($batchq_param_fullfile);
my @jump_q_params = @$first1;
my @nTMTs = @$second1;
if ($nbatch_Q!=$nbatch_internal) {
	die "In ".$ARGV[0].", 'internal_standard_batch' and 'input_n_batch' should have the same batch nos.\n";
}
if (!-f $batchID_fullfile) {
	die "In ".$ARGV[0].", 'path_batch_id' is incorrect (as in 'path_batch_id' there is no -batch-id params file).\n";
}

my @dirs = split('/',$batchID_fullfile); # in_batchidDir, batchID_file
my $in_batchidDir = join('/',@dirs[0..$#dirs-1]); # 'path_batch_id': output path of -batch-id results
my $batchID_file = $dirs[$#dirs];

my ($nbatch_ID,$nbatch_fj_params) = get_nbatch_batchID($batchID_fullfile);
if ($nbatch_ID!=$nbatch_Q) {
	die $ARGV[0]." and ".$batchID_file." should have the same batch nos.\n";
}

# batches
my @bkeys1; # b1, b2, ...
my @bkeys2; # batch1, batch2, ...
for (my $i = 0; $i < $nbatch_Q; $i++) {
	my $ano = $i+1;
	$bkeys1[$i] = "b$ano";
	$bkeys2[$i] = "batch$ano";
}

# quan_dir: $cwd_dir."/".$output_folder_Q."/jump_q_results"
my $qbatch_dir = $cwd_dir."/".$output_folder_Q;
if (!-e $qbatch_dir) {
	system ("mkdir $qbatch_dir");
}
my $quan_dir = $qbatch_dir."/jump_q_results";
if (!-e $quan_dir) {
	system ("mkdir $quan_dir");
}

# batch_paths
my @batch_paths;
for (my $i=0; $i<$#bkeys1+1; $i++){
	$batch_paths[$i] = $quan_dir."/".$bkeys1[$i];
}
make_path(@batch_paths);

# write q params for each batch
my $nchoice4mode3 = 1;# mode3 - peptides: default, from phosphor proteome; else from whole proteome
my $in_ql_fullfile;
if ($input_mode==1) {
	$in_ql_fullfile = $cwd_dir."/ParameterFiles/TMThh/jump_qj_HH_tmt10_human.params";
} elsif ($input_mode==3) {
	if ($nchoice4mode3==1) {# 3: peptides (from phosphor proteome) - default
		$in_ql_fullfile = $cwd_dir."/ParameterFiles/TMThhpho/jump_qj_HH_pho_tmt10_human.params";
	} else {# 3: peptides (from whole proteome)
		$in_ql_fullfile = $cwd_dir."/ParameterFiles/TMThh/jump_qj_HH_tmt10_human.params";
	}
}
my $out_ql_fullfile;
my $maxsize_q = 0;
my $maxsize_l = 0;
my $c_txtsize;
for (my $i = 0; $i < $nbatch_Q; $i++){
	if ($input_mode!=2) {
		$out_ql_fullfile = $batch_paths[$i]."/jump_qj.params";
		write_param2_q($in_ql_fullfile,$out_ql_fullfile,$in_batchidDir,$nTMTs[$i],$bkeys2[$i],$input_mode,$nchoice4mode3,@jump_q_params); # write_param2_q
		$c_txtsize = get_txtsize_from_jump_qj($out_ql_fullfile);
		if ($maxsize_q<$c_txtsize) {
			$maxsize_q = $c_txtsize;
		}
	} else {
		$in_ql_fullfile = $cwd_dir."/ParameterFiles/jump_l.params";
		$out_ql_fullfile = $batch_paths[$i]."/jump_l.params";
		write_param2_l($in_ql_fullfile,$out_ql_fullfile,$in_batchidDir,$bkeys2[$i]); # write_param2_l
		$c_txtsize = get_txtsize_from_jump_lj($out_ql_fullfile);
		if ($maxsize_l<$c_txtsize) {
			$maxsize_l = $c_txtsize;
		}
		
		$in_ql_fullfile = $cwd_dir."/ParameterFiles/TMThhpho/jump_qj_HH_pho_tmt10_human.params";
		$out_ql_fullfile = $batch_paths[$i]."/jump_qj.params";
		write_param2_q($in_ql_fullfile,$out_ql_fullfile,$quan_dir,$nTMTs[$i],$bkeys1[$i],$input_mode,$nchoice4mode3,@jump_q_params); # write_param2_q
		$c_txtsize = get_txtsize_from_jump_qj($out_ql_fullfile);
		if ($maxsize_q<$c_txtsize) {
			$maxsize_q = $c_txtsize;
		}
	}
}

# write_step2_q
my $nMEM1;
my $ql_sh_fullfile;
if ($input_mode==2) {
	$nMEM1 = calculate_mem($maxsize_l,3); # choice: 3,l
	$ql_sh_fullfile = $cwd_dir."/ParameterFiles/Batch/run2_l.sh"; # step2_l (jump -l jump_l.params --queue=standard --mem=20000)
	write_step2_l($ql_sh_fullfile,$quan_dir,$nMEM1,@bkeys1);
	
	# first get mem and queue, and then run step2_l
	submit_jobs($nMEM1,$ql_sh_fullfile);
}

$ql_sh_fullfile = $cwd_dir."/ParameterFiles/Batch/run2_q.sh";
write_step2_q($ql_sh_fullfile,$quan_dir,@bkeys1);

# first get mem and queue, and then run step2_q
$nMEM1 = calculate_mem($maxsize_q,2); # choice: 2,q
submit_jobs($nMEM1,$ql_sh_fullfile);

# write qbatch params (to add or update correct 'input_path_batch')
write_param3_qbatch($batchq_param_fullfile,$input_mode,$nchoice4mode3,@batch_paths);

# write_step3_qbatch
my $qbatch_sh_fullfile = $cwd_dir."/ParameterFiles/Batch/run3_qbatch.sh";
write_step3_qbatch($qbatch_sh_fullfile,$cwd_dir,$Bin,$ARGV[0]);

# first get mem and queue, and then run step3_qbatch
$c_txtsize = get_txtsize_from_qbatch($batchq_param_fullfile);
$nMEM1 = calculate_mem($c_txtsize,5); # choice: 5,qbatch
submit_jobs($nMEM1,$qbatch_sh_fullfile);
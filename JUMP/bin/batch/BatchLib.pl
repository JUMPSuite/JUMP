#!/bin/env perl
use strict;
use warnings;
use autodie; # die if problem reading or writing a file
use Cwd;

sub get_nbatch_batchID {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# nbatch_ID
	my ($first2) = locate_same_head_all("input_path_batch",@mystrs1);
	my @st = @$first2;
	my $nbatch_ID = $#st+1;
	
	# nbatch_fj_params
	my $nbatch_fj_params = 0;
	for (my $i=0; $i<$nbatch_ID; $i++) {
		my $jump_fj_fullfile = get_jump_fj($mystrs1[$st[$i]]);
		if (-e $jump_fj_fullfile) {
			$nbatch_fj_params++;
		}
	}
	
	return $nbatch_ID,$nbatch_fj_params;
}

sub write_step1_fbatch {
	my ($in_sh_fullfile,$in_DataDir,$in_BinDir,$in_param_file) = @_; # input params
	my @mystrs;
	my $lno = 0;
	
	$mystrs[$lno] = "#!/bin/bash\n";
	$lno++;
	$mystrs[$lno] = "cd ".$in_DataDir."\n";
	$lno++;
	$mystrs[$lno] = "module load jump/1.13.003\n"; # include Comet
	$lno++;
	$mystrs[$lno] = "module load perl/5.10.1\n";
	$lno++;
	$mystrs[$lno] = "export PERL5LIB=\"\$PERL5LIB:/home/yli4/lib/perl5/\"\n";
	$lno++;
	$mystrs[$lno] = "perl ".$in_BinDir."/batch/jump_f_batch.pl ".$in_param_file."\n";
	$lno++;
	
	# write
	write_ALL($in_sh_fullfile, @mystrs);
}

sub get_nbatch_batchQ {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# nbatch_Q
	my ($first2) = locate_same_head_all("input_n_batch",@mystrs1);
	my @st = @$first2;
	my $nbatch_Q = $#st+1;
	
	# nbatch_internal
	my ($first3) = locate_same_head_all("internal_standard_batch",@mystrs1);
	my @internal = @$first3;
	my $nbatch_internal = $#internal+1;
	
	# input_mode
	my $st11 = locate_oneline("input_mode",@mystrs1);
	my $value1 = get_value($mystrs1[$st11]);
	my $input_mode = int($value1);
	
	# batchID_fullfile
	my $st12 = locate_oneline("path_batch_id",@mystrs1);
	my $path_batch_id = get_value($mystrs1[$st12]);
	my $batchID_fullfile = get_theonlyone_params($path_batch_id);
	
	# output_folder_Q: output folder suffix name; prefix always 'batch'
	my $st13 = locate_oneline("output_folder",@mystrs1);
	my $value2 = get_value($mystrs1[$st13]);
	my $output_folder_Q = "batch_".$value2;
	
	# jump_q_params
	my $st14 = locate_oneline("ppi_filter",@mystrs1);
	my $st15 = locate_oneline("impurity_correction",@mystrs1);
	my $st16 = locate_oneline("loading_bias_correction",@mystrs1);
	my $st17 = locate_oneline("interference_removal",@mystrs1);
	my @jump_q_params;
	my $jno = 0;
	$jump_q_params[$jno] = $mystrs1[$st14];
	$jno++;
	$jump_q_params[$jno] = $mystrs1[$st15];
	$jno++;
	$jump_q_params[$jno] = $mystrs1[$st16];
	$jno++;
	$jump_q_params[$jno] = $mystrs1[$st17];
	$jno++;
	
	# nTMTs
	my @nTMTs;
	for (my $i=0; $i<$nbatch_Q; $i++) {
		my $value2 = get_value($mystrs1[$st[$i]]);
		$nTMTs[$i] = int($value2);
	}
	
	return $nbatch_Q,$nbatch_internal,$input_mode,$batchID_fullfile,$output_folder_Q,\@jump_q_params,\@nTMTs;
}

sub get_value {
	my ($in_str) = @_; # 1st input param
	if (index($in_str, '#')!=-1) { # it has '#' (*preferred, if no '#', the val is with '\n', which will be transferred to later replace use*)
		my ($L,$L0) = split(/#/,$in_str);
		$in_str = $L;
	}
	my ($mykey,$myvalue) = split(/=/,$in_str);
	$myvalue =~ s/\s+//g; # remove [\ \t\r\n\f]
	
	return $myvalue;
}

sub get_theonlyone_params {
	my ($ID_dir) = @_;
	
	my $jump_1params_fullfile = '';
	opendir(my $dh, $ID_dir) or die $!;
	while (my $file = readdir($dh)) {
		my $myfullfile = $ID_dir."/".$file;
		if (-f $myfullfile && index($file, '.params')!=-1) {
			$jump_1params_fullfile = $myfullfile;
			last;
		}
	}
	closedir($dh);
	
	return $jump_1params_fullfile;
}

sub write_param2_q {
	my ($in_q_fullfile,$out_q_fullfile,$in_batchidDir,$in_nTMT,$in_batch1,$input_mode,$nchoice4mode3,@jump_q_params) = @_; # input params
	my ($first1) = read_ALL($in_q_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# locate and replace
	my $st11 = locate_oneline("idtxt",@mystrs1);
	if ($input_mode==1) {
		$mystrs1[$st11] = "idtxt = ".$in_batchidDir."/".$in_batch1."/sum_out/ID.txt\n";
	} elsif ($input_mode==3) {
		if ($nchoice4mode3==1) {
			$mystrs1[$st11] = "idtxt = ".$in_batchidDir."/".$in_batch1."/sum_out_mod/IDmod.txt\n";
		} else {
			$mystrs1[$st11] = "idtxt = ".$in_batchidDir."/".$in_batch1."/sum_out/ID.txt\n";
		}
	} else {
		$mystrs1[$st11] = "idtxt = ".$in_batchidDir."/".$in_batch1."/loc_out_mod/IDmod.txt\n";
	}
	
	# locate and replace
	my $st12 = locate_oneline("save_dir",@mystrs1);
	if ($input_mode==1) {
		$mystrs1[$st12] = "save_dir = HH_tmt							# name of the directory for JUMPq results (prefix \"quan-\" will be added)\n";
	} else {
		if ($input_mode==3 and $nchoice4mode3!=1) {
			$mystrs1[$st12] = "save_dir = HH_tmt							# name of the directory for JUMPq results (prefix \"quan-\" will be added)\n";
		} else {
			$mystrs1[$st12] = "save_dir = HH_tmt_mod							# name of the directory for JUMPq results (prefix \"quan-\" will be added)\n";
		}
	}
	
	# create_TMT4q
	my ($impurity_matrix, $tmt_reporters_used, $comparison_groups_twoGroups, @sample_labels) = create_TMT4q($in_nTMT);
	
	# locate and replace
	my $st13 = locate_oneline("impurity_matrix",@mystrs1);
	$mystrs1[$st13] = $impurity_matrix;
	
	# locate and replace
	my $st14 = locate_oneline("tmt_reporters_used",@mystrs1);
	$mystrs1[$st14] = $tmt_reporters_used;
	
	# comparison_analysis
	my $st015 = locate_oneline("comparison_analysis",@mystrs1);
	my $value1 = get_value($mystrs1[$st015]);
	my $comparison_analysis = int($value1);
	
	# locate and replace
	my $st15 = locate_oneline("comparison_groups_twoGroups",@mystrs1);
	if ($comparison_analysis==0) {
		$mystrs1[$st15] = $comparison_groups_twoGroups;
	}
	
	# the values from qbatch will overwrite the default values
	my $st16 = locate_oneline("ppi_filter",@mystrs1);
	my $st17 = locate_oneline("impurity_correction",@mystrs1);
	my $st18 = locate_oneline("loading_bias_correction",@mystrs1);
	my $st19 = locate_oneline("interference_removal",@mystrs1);
	$mystrs1[$st16] = $jump_q_params[0];
	$mystrs1[$st17] = $jump_q_params[1];
	$mystrs1[$st18] = $jump_q_params[2];
	$mystrs1[$st19] = $jump_q_params[3];
	
	# locate
	my ($st1,$st2) = locate_same_head('sig',@mystrs1);
	
	# merge
	my ($first2) = merge_one_part($st1,$st2,\@mystrs1,\@sample_labels); # (*passing multi arrays is the same way as returning*)
	@mystrs1 = @$first2;
	
	# write
	write_ALL($out_q_fullfile, @mystrs1);
}

sub write_param2_l {
	my ($in_l_fullfile,$out_l_fullfile,$in_batchidDir,$in_batch1) = @_; # input params
	my ($first1) = read_ALL($in_l_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# locate and replace
	my $st11 = locate_oneline("IDmod",@mystrs1);
	$mystrs1[$st11] = "IDmod = ".$in_batchidDir."/".$in_batch1."/sum_out_mod/IDmod.txt\n";
	
	# write
	write_ALL($out_l_fullfile, @mystrs1);
}

sub write_param3_qbatch {
	my ($in_qbatch_fullfile,$input_mode,$nchoice4mode3,@batch_paths) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_qbatch_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# input_path_batch
	my @input_path_batch;
	my $nbatch = $#batch_paths+1;
	for (my $i=0; $i<$nbatch; $i++) {
		my $ino = $i+1;
		if ($input_mode==1) {
			$input_path_batch[$i] = "input_path_batch".$ino." = ".$batch_paths[$i]."/quan_HH_tmt/publications\n";
		} else {
			if ($input_mode==3 and $nchoice4mode3!=1) {
				$input_path_batch[$i] = "input_path_batch".$ino." = ".$batch_paths[$i]."/quan_HH_tmt/publications\n";
			} else {
				$input_path_batch[$i] = "input_path_batch".$ino." = ".$batch_paths[$i]."/quan_HH_tmt_mod/publications\n";
			}
		}
	}
	
	# locate
	my ($st1,$st2) = locate_same_head('input_path_batch',@mystrs1);
	
	# add or update correct 'input_path_batch'
	if ($st1==-1) { # no 'input_path_batch' in jump_batchQ.params
		my $nparams_lines = $#mystrs1+1;
		$mystrs1[$nparams_lines] = "\n";
		$mystrs1[$nparams_lines+1] = "# Inputs: absolute path of publication tables from JUMP -q results\n";
		for (my $i=0; $i<$nbatch; $i++) {
			$mystrs1[$i+$nparams_lines+2] = $input_path_batch[$i];
		}
	} else { # it already has 'input_path_batch' in jump_batchQ.params
		my ($first2) = merge_one_part($st1,$st2,\@mystrs1,\@input_path_batch);
		@mystrs1 = @$first2;
	}
	
	# write
	write_ALL($in_qbatch_fullfile, @mystrs1);
}

sub write_step2_q {
	my ($in_sh_fullfile,$in_QuanDir,@in_bs) = @_; # input params
	my @mystrs;
	my $lno = 0;
	
	$mystrs[$lno] = "#!/bin/bash\n";
	$lno++;
	$mystrs[$lno] = "cd ".$in_QuanDir."\n";
	$lno++;
	$mystrs[$lno] = "module load jump/1.13.003\n";
	$lno++;
	for (my $i=0; $i<$#in_bs+1; $i++) {
		$mystrs[$lno] = "cd ".$in_bs[$i]."\n";
		$lno++;
		$mystrs[$lno] = "jump -q jump_qj.params -dispatch=localhost\n";
		$lno++;
		$mystrs[$lno] = "cd ..\n";
		$lno++;
	}
	
	# write
	write_ALL($in_sh_fullfile, @mystrs);
}

sub write_step2_l {
	my ($in_sh_fullfile,$in_QuanDir,$nMEM1,@in_bs) = @_; # input params
	my $default_queue = "standard";
	if ($nMEM1>=200) {
		$default_queue = "large_mem";
	}
	my $nMEM2 = $nMEM1*1000;
	my @mystrs;
	my $lno = 0;
	
	$mystrs[$lno] = "#!/bin/bash\n";
	$lno++;
	$mystrs[$lno] = "cd ".$in_QuanDir."\n";
	$lno++;
	$mystrs[$lno] = "module load jump/1.13.003\n";
	$lno++;
	for (my $i=0; $i<$#in_bs+1; $i++) {
		$mystrs[$lno] = "cd ".$in_bs[$i]."\n";
		$lno++;
		$mystrs[$lno] = "jump -l jump_l.params --queue=".$default_queue." --mem=".$nMEM2."\n";
		$lno++;
		$mystrs[$lno] = "cd ..\n";
		$lno++;
	}
	
	# write
	write_ALL($in_sh_fullfile, @mystrs);
}

sub write_step3_qbatch {
	my ($in_sh_fullfile,$in_DataDir,$in_BinDir,$in_param_file) = @_; # input params
	my @mystrs;
	my $lno = 0;
	
	$mystrs[$lno] = "#!/bin/bash\n";
	$lno++;
	$mystrs[$lno] = "cd ".$in_DataDir."\n";
	$lno++;
	$mystrs[$lno] = "module load jump/1.13.003\n"; # include Comet
	$lno++;
	$mystrs[$lno] = "module load perl/5.10.1\n";
	$lno++;
	$mystrs[$lno] = "export PERL5LIB=\"\$PERL5LIB:/home/yli4/lib/perl5/\"\n";
	$lno++;
	$mystrs[$lno] = "perl ".$in_BinDir."/batch/jump_batch_v1.6.pl ".$in_param_file."\n";
	$lno++;
	
	# write
	write_ALL($in_sh_fullfile, @mystrs);
}

sub submit_jobs {
	my ($nMEM1,$sh_fullfile) = @_; # 1st input param
	my $default_queue = "standard";
	if ($nMEM1>=200) {
		$default_queue = "large_mem";
	}
	my $nMEM2 = $nMEM1*1000;
	my $hint1 = "Applying ".$nMEM1." GB RAM in queue <".$default_queue."> (please be patient)\n";
	my $cmd1 = "bsub -R \"rusage[mem=".$nMEM2."]\" -q ".$default_queue." -P Proteomics -Is bash ".$sh_fullfile;
	print $hint1;
	system($cmd1);
}

## get_pepXMLsize_from_jump_fj
# $pepXMLsize = 0;
# $engine = get_engine_from_jump_fj($jump_fj_fullfile);
# $pepXMLsize += get_pepXMLsize_from_jump_fj($jump_fj_fullfile,$engine);

sub get_txtsize_from_jump_qj {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# idtxt
	my $st11 = locate_oneline("idtxt",@mystrs1);
	my $idtxt_fullfile = get_value($mystrs1[$st11]);
	my $txtsize = 0;
	if (-f $idtxt_fullfile) {
		$txtsize += (stat $idtxt_fullfile)[7]/(1024*1024);
	}
	
	return $txtsize;
}

sub get_txtsize_from_jump_lj {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# IDmod
	my $st11 = locate_oneline("IDmod",@mystrs1);
	my $idtxt_fullfile = get_value($mystrs1[$st11]);
	my $txtsize = 0;
	if (-f $idtxt_fullfile) {
		$txtsize += (stat $idtxt_fullfile)[7]/(1024*1024);
	}
	
	return $txtsize;
}

sub get_nfraction_from_fbatch {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# nbatch_ID
	my ($first2) = locate_same_head_all("input_path_batch",@mystrs1);
	my @st = @$first2;
	my $nbatch_ID = $#st+1;
	
	# nfraction, pepXMLsize
	my $nfraction = 0;
	for (my $i=0; $i<$nbatch_ID; $i++) {
		my $jump_fj_fullfile = get_jump_fj($mystrs1[$st[$i]]);
		if (-e $jump_fj_fullfile) {
			$nfraction += get_nfraction_from_jump_fj($jump_fj_fullfile);
		}
	}
	
	return $nfraction;
}

sub get_pepXMLsize_from_fbatch {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# nbatch_ID
	my ($first2) = locate_same_head_all("input_path_batch",@mystrs1);
	my @st = @$first2;
	my $nbatch_ID = $#st+1;
	
	# nfraction, pepXMLsize
	my $pepXMLsize = 0;
	for (my $i=0; $i<$nbatch_ID; $i++) {
		my $jump_fj_fullfile = get_jump_fj($mystrs1[$st[$i]]);
		my $engine = get_engine_from_jump_fj($jump_fj_fullfile);
		if (-e $jump_fj_fullfile) {
			$pepXMLsize += get_pepXMLsize_from_jump_fj($jump_fj_fullfile,$engine);
		}
	}
	
	return $pepXMLsize;
}

sub get_txtsize_from_qbatch {
	my ($in_param_fullfile) = @_; # input params (*when mix of scalar and array, array is the last*)
	my $nlevel = 1;
	#my ($in_param_fullfile,$nlevel) = @_; # input params (*when mix of scalar and array, array is the last*)
	my ($first1) = read_ALL($in_param_fullfile); # read the entire param file
	my @mystrs1 = @$first1;
	
	# nbatch_ID
	my ($first2) = locate_same_head_all("input_path_batch",@mystrs1);
	my @st = @$first2;
	my $nbatch_ID = $#st+1;
	
	# filename1, filename2
	my $filename1;
	my $filename2;
	if ($nlevel==1) {
		$filename1 = "id_all_pep_quan.txt";
		$filename2 = "id_uni_pep_quan.txt";
	} else {
		$filename1 = "id_all_prot_quan.txt";
		$filename2 = "id_uni_prot_quan.txt";
	}
	
	# txtsize
	my $txtsize = 0;
	for (my $i=0; $i<$nbatch_ID; $i++) {
		my $input_path_batch = get_value($mystrs1[$st[$i]]);
		my $fullname1 = $input_path_batch."/".$filename1;
		my $fullname2 = $input_path_batch."/".$filename2;
		if (-f $fullname1) {
			$txtsize += (stat $fullname1)[7]/(1024*1024);
		}
		if (-f $fullname2) {
			$txtsize += (stat $fullname2)[7]/(1024*1024);
		}
	}
	
	return $txtsize;
}

sub calculate_mem {
	my ($inputsize,$choice) = @_; # input params
	# choice: 1,f; 2,q; 3,l; 4,fbatch; 5,qbatch
	
	my $slope;
	my $intercept;
	if ($choice==1) {# choice: 1,f
		$slope = 8;
		$intercept = 550;
	} elsif ($choice==2) {# choice: 2,q
		$slope = 62.4775;
		$intercept = 378.6;
	} elsif ($choice==3) {# choice: 3,l (the same as 2,q)
		$slope = 62.4775;
		$intercept = 378.6;
	} elsif ($choice==4) {# choice: 4,fbatch
		$slope = 8.6277;
		$intercept = 2798.3;
	} else {# choice: 5,qbatch
		$slope = 8.1035;
		$intercept = 41225.6;
	}
	my $nMEM0 = int( ($slope*$inputsize+$intercept)/1000 + 0.5 ); # from inputsize to size of memory (GB)
	my $tol = 10;
	my $nMEM1 = $nMEM0+$tol; # add 10 GB;
	
	return $nMEM1;
}






sub read_ALL {
	my ($in_fullfile) = @_; # 1st input param
	my @mystrs; # arrays to return
	my $kno = 0; # index of arrays
	open my $fh, '<:encoding(utf8)', $in_fullfile;
	while (my $line = <$fh>) { # get a line
		$mystrs[$kno] = $line;
		$kno++;
	}
	close($fh);
	return \@mystrs;
}

sub write_ALL {
	my ($in_fullfile, @mystrs) = @_; # 1st input param
	open my $fh, '>:encoding(utf8)', $in_fullfile;
	foreach my $line ( @mystrs ) { # add the line to the file
		print {$fh} $line;
	}
	close($fh);
}

sub locate_same_head_all {
	my ($head,@mystrs) = @_; # input params
	my @st;
	my $sno = 0;
	for (my $i = 0; $i < $#mystrs+1; $i++) { # locate start
		my $c_str = $mystrs[$i];
		$c_str =~ s/\s+//g; # remove [\ \t\r\n\f]
		if (index($c_str, $head)==0) {
			$st[$sno] = $i;
			$sno++;
		}
	}
	return \@st;
}

sub locate_same_head {
	my ($head,@mystrs) = @_; # input params
	my $st1 = -1;
	my $st2 = -1;
	for (my $i = 0; $i < $#mystrs+1; $i++) { # locate start
		my $c_str = $mystrs[$i];
		$c_str =~ s/\s+//g; # remove [\ \t\r\n\f]
		if (index($c_str, $head)==0) {
			$st1 = $i;
			last;
		}
	}
	$st2 = $st1;
	for (my $i = $st1+1; $i < $#mystrs+1; $i++) { # locate terminate
		my $c_str = $mystrs[$i];
		$c_str =~ s/\s+//g; # remove [\ \t\r\n\f]
		if (index($c_str, $head)==0) {
			$st2 = $i;
		}
	}
	return $st1, $st2;
}

sub get_jump_fj {
	my ($in_str) = @_;
	my ($mykey,$ID_fullfile) = split(/=/,$in_str);
	$ID_fullfile =~ s/\s+//g; # remove [\ \t\r\n\f]
	my @dirs = split('/',$ID_fullfile);
	my $ID_dir = join('/',@dirs[0..$#dirs-1]);
	
	my $jump_fj_fullfile = '';
	opendir(my $dh, $ID_dir) or die $!;
	while (my $file = readdir($dh)) {
		my $myfullfile = $ID_dir."/".$file;
		if (-f $myfullfile && index($file, '.params')!=-1) {
			$jump_fj_fullfile = $myfullfile;
			last;
		}
	}
	closedir($dh);
	
	return $jump_fj_fullfile;
}

sub get_engine_from_jump_fj {
	my ($in_fullfile) = @_; # input params
	my ($first1) = read_ALL($in_fullfile); # read the entire param file
	my @mystrs = @$first1;
	
	# engine
	my $key = "search_engine = comet";
	$key =~ s/\s+//g; # remove [\ \t\r\n\f]
	my $engine = 1; # 1: JUMP; 2: Comet
	if ( locate_oneline($key,@mystrs)!=-1 ) {
		$engine = 2;
	}
	
	return $engine;
}

sub get_nfraction_from_jump_fj {
	my ($in_fullfile) = @_; # input params
	my ($first1) = read_ALL($in_fullfile); # read the entire param file
	my @mystrs = @$first1;
	
	# locate
	my ($st1,$st2) = locate11(@mystrs);
	
	# count
	my $nfraction = 0;
	for (my $i = $st1; $i <= $st2; $i++) {
		if ( index($mystrs[$i], '/')!=-1 ) {
			my $ID_fullpath = $mystrs[$i];
			$ID_fullpath =~ s/\s+//g; # remove [\ \t\r\n\f]
			if (index($ID_fullpath, ':')!=-1) { # has HH_tmt: in the beginning
				my ($head,$path) = split(':',$ID_fullpath);
				$ID_fullpath = $path;
			}
			if (index($ID_fullpath, '#')==0) { # skip the one excluded
				next;
			}
			
			$nfraction++;
		}
	}
	return $nfraction;
}

sub get_pepXMLsize_from_jump_fj {
	my ($in_fullfile,$in_engine) = @_; # input params
	my ($first1) = read_ALL($in_fullfile); # read the entire param file
	my @mystrs = @$first1;
	
	# locate
	my ($st1,$st2) = locate11(@mystrs);
	
	# count
	my @pepxml = (".pepXML",".pep.xml");# 1: JUMP; 2: Comet
	my $pepXMLsize = 0;
	for (my $i = $st1; $i <= $st2; $i++) {
		if ( index($mystrs[$i], '/')!=-1 ) {
			my $ID_fullpath = $mystrs[$i];
			$ID_fullpath =~ s/\s+//g; # remove [\ \t\r\n\f]
			if (index($ID_fullpath, ':')!=-1) { # has HH_tmt: in the beginning
				my ($head,$path) = split(':',$ID_fullpath);
				$ID_fullpath = $path;
			}
			if (index($ID_fullpath, '#')==0) { # skip the one excluded
				next;
			}
			
			# my @dirs = split('/',$ID_fullpath); # /NY198case_b01_f01/NY198case_b01_f01.1 (JUMP) or /NCI-11plex-1-F1-f10268/NCI-11plex-1-F1-f10268.1 (Comet)
			# my $ID_dir = join('/',@dirs[0..$#dirs-1]); # /NY198case_b01_f01 or /NCI-11plex-1-F1-f10268
			# # if /NY198case_b01_f01/NY198case_b01_f01.1/ (JUMP) or /NCI-11plex-1-F1-f10268/NCI-11plex-1-F1-f10268.1/ (Comet)
			# if (index($dirs[$#dirs-1], $dirs[$#dirs-2])==0) {
				# $ID_dir = join('/',@dirs[0..$#dirs-2]);
			# }
			
			# my $pepxml_fullfile = '';
			# opendir(my $dh, $ID_dir) or die $!;
			# while (my $file = readdir($dh)) {
				# my $myfullfile = $ID_dir."/".$file;
				# if (-f $myfullfile && index($file, $pepxml[$in_engine-1])!=-1) {
					# $pepxml_fullfile = $myfullfile;
					# last;
				# }
			# }
			# closedir($dh);
			
			my @dirs = split('/',$ID_fullpath); # /NY198case_b01_f01/NY198case_b01_f01.1 (JUMP) or /NCI-11plex-1-F1-f10268/NCI-11plex-1-F1-f10268.1 (Comet)
			# if /NY198case_b01_f01/NY198case_b01_f01.1/ (JUMP) or /NCI-11plex-1-F1-f10268/NCI-11plex-1-F1-f10268.1/ (Comet)
			$ID_fullpath = join('/',@dirs[0..$#dirs]);# it will get rid of the last '/' if it is there
			my $pepxml_fullfile = $ID_fullpath.$pepxml[$in_engine-1];
			
			if (-f $pepxml_fullfile) {
				$pepXMLsize += (stat $pepxml_fullfile)[7]/(1024*1024);
			}
		}
	}
	return $pepXMLsize;
}

sub locate_oneline {
	my ($key,@mystrs) = @_; # input params
	my $st1 = -1;
	for (my $i = 0; $i < $#mystrs+1; $i++) { # locate start
		my $c_str = $mystrs[$i];
		$c_str =~ s/\s+//g; # remove [\ \t\r\n\f]
		if (index($c_str, $key)==0) {
			$st1 = $i;
			last;
		}
	}
	return $st1;
}

sub locate11 {
	my (@mystrs) = @_; # input params
	my $st1 = -1;
	my $st2 = -1;
	for (my $i = 0; $i < $#mystrs+1; $i++) { # locate start
		my $c_str = $mystrs[$i];
		$c_str =~ s/\s+//g; # remove [\ \t\r\n\f]
		if ( index($c_str, ':/')!=-1 && !(substr($c_str,0,1) eq "#") ) {
			$st1 = $i;
			last;
		}
	}
	$st2 = $st1;
	for (my $i = $st1+1; $i < $#mystrs+1; $i++) { # locate terminate
		my $c_str = $mystrs[$i];
		$c_str =~ s/\s+//g; # remove [\ \t\r\n\f]
		if (index($c_str, '/')==0) {
			$st2 = $i;
		}
	}
	return $st1, $st2;
}

sub create_TMT4q {
	my ($nTMT) = @_; # input params
	my $impurity_matrix;
	my $tmt_reporters_used;
	my $comparison_groups_twoGroups;
	my @sample_labels;
	my $ini_file;
	if ($nTMT==6 || $nTMT==8 || $nTMT==10 || $nTMT==11 || $nTMT==16) {
		$ini_file = "TMT".$nTMT.".ini";
	}
	$impurity_matrix = "impurity_matrix = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/JUMPq/".$ini_file."	# impurity table for correction\n";
	
	if ($nTMT==6) {
		$tmt_reporters_used = "tmt_reporters_used = sig126; sig127; sig128; sig129; sig130; sig131\n";
		$comparison_groups_twoGroups = "comparison_groups_twoGroups = S1, S2, S3 : S4, S5, S6\n";
		$sample_labels[0] = "sig126 = S1\n";
		$sample_labels[1] = "sig127 = S2\n";
		$sample_labels[2] = "sig128 = S3\n";
		$sample_labels[3] = "sig129 = S4\n";
		$sample_labels[4] = "sig130 = S5\n";
		$sample_labels[5] = "sig131 = S6\n";
	} elsif ($nTMT==8) {
		$tmt_reporters_used = "tmt_reporters_used = sig126; sig127N; sig127C; sig128; sig129N; sig129C; sig130; sig131\n";
		$comparison_groups_twoGroups = "comparison_groups_twoGroups = S1, S2, S3, S4 : S5, S6, S7, S8\n";
		$sample_labels[0] = "sig126 = S1\n";
		$sample_labels[1] = "sig127N = S2\n";
		$sample_labels[2] = "sig127C = S3\n";
		$sample_labels[3] = "sig128 = S4\n";
		$sample_labels[4] = "sig129N = S5\n";
		$sample_labels[5] = "sig129C = S6\n";
		$sample_labels[6] = "sig130 = S7\n";
		$sample_labels[7] = "sig131 = S8\n";
	} elsif ($nTMT==10) {
		$tmt_reporters_used = "tmt_reporters_used = sig126; sig127N; sig127C; sig128N; sig128C; sig129N; sig129C; sig130N; sig130C; sig131\n";
		$comparison_groups_twoGroups = "comparison_groups_twoGroups = S1, S2, S3, S4, S5 : S6, S7, S8, S9, S10\n";
		$sample_labels[0] = "sig126 = S1\n";
		$sample_labels[1] = "sig127N = S2\n";
		$sample_labels[2] = "sig127C = S3\n";
		$sample_labels[3] = "sig128N = S4\n";
		$sample_labels[4] = "sig128C = S5\n";
		$sample_labels[5] = "sig129N = S6\n";
		$sample_labels[6] = "sig129C = S7\n";
		$sample_labels[7] = "sig130N = S8\n";
		$sample_labels[8] = "sig130C = S9\n";
		$sample_labels[9] = "sig131 = S10\n";
	} elsif ($nTMT==11) {
		$tmt_reporters_used = "tmt_reporters_used = sig126; sig127N; sig127C; sig128N; sig128C; sig129N; sig129C; sig130N; sig130C; sig131N; sig131C\n";
		$comparison_groups_twoGroups = "comparison_groups_twoGroups = S1, S2, S3, S4, S5 : S6, S7, S8, S9, S10, S11\n";
		$sample_labels[0] = "sig126 = S1\n";
		$sample_labels[1] = "sig127N = S2\n";
		$sample_labels[2] = "sig127C = S3\n";
		$sample_labels[3] = "sig128N = S4\n";
		$sample_labels[4] = "sig128C = S5\n";
		$sample_labels[5] = "sig129N = S6\n";
		$sample_labels[6] = "sig129C = S7\n";
		$sample_labels[7] = "sig130N = S8\n";
		$sample_labels[8] = "sig130C = S9\n";
		$sample_labels[9] = "sig131N = S10\n";
		$sample_labels[10] = "sig131C = S11\n";
	} elsif ($nTMT==16) {
		$tmt_reporters_used = "tmt_reporters_used = sig126; sig127N; sig127C; sig128N; sig128C; sig129N; sig129C; sig130N; sig130C; sig131N; sig131C; sig132N; sig132C; sig133N; sig133C; sig134N\n";
		$comparison_groups_twoGroups = "comparison_groups_twoGroups = S1, S2, S3, S4, S5, S6, S7, S8 : S9, S10, S11, S12, S13, S14, S15, S16\n";
		$sample_labels[0] = "sig126 = S1\n";
		$sample_labels[1] = "sig127N = S2\n";
		$sample_labels[2] = "sig127C = S3\n";
		$sample_labels[3] = "sig128N = S4\n";
		$sample_labels[4] = "sig128C = S5\n";
		$sample_labels[5] = "sig129N = S6\n";
		$sample_labels[6] = "sig129C = S7\n";
		$sample_labels[7] = "sig130N = S8\n";
		$sample_labels[8] = "sig130C = S9\n";
		$sample_labels[9] = "sig131N = S10\n";
		$sample_labels[10] = "sig131C = S11\n";
		$sample_labels[11] = "sig132N = S12\n";
		$sample_labels[12] = "sig132C = S13\n";
		$sample_labels[13] = "sig133N = S14\n";
		$sample_labels[14] = "sig133C = S15\n";
		$sample_labels[15] = "sig134N = S16\n";
	}
	return $impurity_matrix, $tmt_reporters_used, $comparison_groups_twoGroups, @sample_labels;
}

sub merge_one_part {
	my ($st1,$st2,$first,$second) = @_; # input params
	my @mystrs1 = @$first;
	my @mystrs2 = @$second;
	my @mystrs;
	my $mno = 0;
	for (my $i1=0; $i1<=$st1-1; $i1++){
		$mystrs[$mno] = $mystrs1[$i1];
		$mno++;
	}
	for (my $i2=0; $i2<$#mystrs2+1; $i2++){
		$mystrs[$mno] = $mystrs2[$i2];
		$mno++;
	}
	
	for (my $i3=$st2+1; $i3<$#mystrs1+1; $i3++){
		$mystrs[$mno] = $mystrs1[$i3];
		$mno++;
	}
	return \@mystrs;
}

1;
# Note that you need the 1; at the end of the file.
# This is because Perl needs the last expression in the file to return a true value.
#!/usr/bin/perl 

use strict;
use warnings;
use File::Basename;
use pepXML_parser;
use Clone qw(clone);
use Storable;
use Cwd;

use IDtxt_parser;
my $ip=IDtxt_parser->new;

if (scalar(@ARGV)!=1)
{
        die "Usage: perl jump_f_batch.pl jump_fbatch.params\n";
}

# part 1: parse .params
# parse .params
my (%parahash);
$ip->parse_params($ARGV[0],\%parahash);
$parahash{score_cutoff_quantile}=$parahash{score_cutoff_quantile}/100;

# mkdir
if (!(-e "$parahash{output_folder}")) {
	system(qq(mkdir -p $parahash{output_folder}/intermediate));
}
system(qq(cp $ARGV[0] $parahash{output_folder})); # cp jump_fbatch.params to output folder

my $log_file="$parahash{output_folder}/jump_f_batch.log";
open(LOG,">$log_file");
close LOG;

# %parahash{batches}{batch1/2}=path1/2
foreach my $k (keys %parahash) {
	next unless ($k =~ m/input_path_batch(\d+)/);
	my $b="batch$1";
	$parahash{batches}{$b}=$parahash{$k};
}

# part 2: parse peptide and PSM information for each batch
my (%pep_acp, %pep_fine, %pep2pro, %prohash, %scoreCutoff, %out2pep, %pepXML_outhash, %pepXML_parahash,%pepxmlPath, %jumpf_outfiles, $searchEngine);
foreach my $b (sort keys %{$parahash{batches}}) {
	printLog("  \nAnalyzing $b:\n");

	# section a)
	# parse jump_f.params from original jump -f output
	my $jump_f_params=dirname($parahash{batches}{$b}); # /research/rgs01/project_space/penggrp/BrainTumor_JP/penggrp/test_batchf/sum_Batch1/IDwDecoy.txt
	$jump_f_params.='/jump_f.params';
	printLog("  Parsing jump_f.params ($jump_f_params)\n");
	my %jump_f_parahash;
	$ip->parse_params($jump_f_params,\%jump_f_parahash);
	$searchEngine='jump';
	if (defined($jump_f_parahash{search_engine})) {
		$searchEngine=$jump_f_parahash{search_engine};
		printLog("  Search engine: $searchEngine\n");
	}
	# parse ID.txt
	# use IDtxt_parser->parse_IDtxt()
	printLog("  Parsing ID.txt ($parahash{batches}{$b})\n");
	my ($scanhash,$peptidehash)=$ip->parse_IDtxt($parahash{batches}{$b},1);

	# record -f param path (use f.params as a template for re-running -f in the end)
	#$parahash{paths}{f_params}=dirname($parahash{batches}{$b});
	#$parahash{paths}{f_params}.='/jump';

	# get .pepXML path: 
	# infor from ID.txt: /home/bbai/ADproject/mouse/pho/p39/p39.1/p39.83007.1.3.spout
	# pepXML path: /home/bbai/ADproject/mouse/pho/p39/p39.1.pepXML
	# record .pepXML path: 
	# %pepxmlPath{batch1}{./home/bbai/ADproject/mouse/pho/p39/p39.1.pepXML.}
	#my %pepxmlPath;
	foreach my $outfile (keys %$scanhash) {
		my $pth=dirname($outfile);
		$pth =~ s/\/$//;
		my $pepxml; # pepxml path
		if ($searchEngine eq 'comet') {
			 $pepxml="$pth\.pep.xml";
		} else {
			 $pepxml="$pth\.pepXML";
		}
		$pepxmlPath{$b}{$pepxml}=1;

		my $out=basename($outfile);
		$out =~ s/\.spout$//;
		$out =~ s/\.out$//;
		$jumpf_outfiles{$out}=1;
	}
	# check path result (for debug)
	#foreach my $pepxml (keys %pepxmlPath) { print "$pepxml\n"; }

	# protein FDR
	# parse id_uni_pep: %pep2pro{$pep}=$protein
	# build %prohash{$pro}{$batch}{$pep}{TD}
	#my $id_uni_pep=dirname($parahash{batches}{$b}); # -f folder path
	#$id_uni_pep.='/publications/id_uni_pep.txt';
	#printLog("  Parsing id_uni_pep.txt ($id_uni_pep)\n");
	#proInfor($id_uni_pep,\%pep2pro,\%prohash,$b);

	my $fPephashPath=dirname($parahash{batches}{$b}); # -f folder path
	$fPephashPath.='/misc/peptidehash.hash';
	printLog("  Parsing -f peptidehash ($fPephashPath)\n");
	proInfor($fPephashPath,\%pep2pro,\%prohash,$b);

	# add this batch specific hash to %pep_acp
	# simplify/format %peptidehash: new function: add_idtxtPephash2pephash
	# new structure: %pep_acp{$pep}{batch1}{$outfile}{Jscore/nbPep/TD/scan/fraction}
	add_idtxtPephash2pephash($b, $peptidehash, $scanhash, \%pep_acp);

	# section b)
	# parse .pepXML
	my $pepxps=pepXML_parser->new;
	my %whlBatchOut;
	foreach my $pepxml (keys %{$pepxmlPath{$b}}) {
		# use pepXML_parser->pepXML2Hashes()
		# get %outInforHash
		printLog("  Parsing .pepXML files ($pepxml) and filter by Jscore of $parahash{min_Jscore} as initial threshold\n");
		#printLog("  Using minimum Jscore of $parahash{min_Jscore} as initial threshold\n");
		my (%paraHash, %outInforHash);
		$pepxps->pepXML2Hashes(\%paraHash, \%outInforHash, $pepxml);

		# add this batch specific hash to %pep_fine
		# convert outfile hash to peptide hash: write a new function: pepxmlOuthash2pephash
		# Jscore cutoff: as determined in parameter
		pepxmlOuthash2pephash($b, \%paraHash, \%outInforHash, \%pep_fine, $parahash{min_Jscore},$parahash{multiHit_max_dJn}, \%out2pep, $searchEngine);

		# build %pepXML_outhash: add outfiles that pass Jscore cutoff
		#buildOutHash(\%whlBatchOut, \%paraHash, \%outInforHash, $parahash{min_Jscore});
		buildOutHash(\%whlBatchOut, \%paraHash, \%outInforHash, $parahash{min_Jscore},\%jumpf_outfiles,$searchEngine);

		# save %paraHash
		$pepXML_parahash{$b}=clone(\%paraHash);

		# record one pepXML path (use to cp one -s output folder)
		$parahash{paths}{pepXML}=$pepxml;

		# delete outhash
		undef %outInforHash;
	}

	$pepXML_outhash{$b}=clone(\%whlBatchOut);

}

# section c): estimate FDR for overlapped/batch-specific PSMs/peptides
#my $nBatch=scalar(keys %{$parahash{batches}});
my (%summaryTable);
$summaryTable{original}{initial}=1;
$summaryTable{batch_improved}{initial}=1;
printLog("\n--------------------------------------- Before PSM rescue ---------------------------------------\n");
printFDR(\%pep_acp, \%prohash, $parahash{batches},1,1,1,1,$summaryTable{original});

# section d): calculate Jscore thresholds for each PSM group defined by charge states and peptide length
groupScore(\%pep_acp,\%scoreCutoff,$parahash{score_cutoff_quantile});
printGroupScoreCuts(\%scoreCutoff,"$parahash{output_folder}/intermediate/group_score_cutoffs.txt");

# part 3: rescue the missing peptides 
#printLog("\n--------------------------------------- Rescuing PSMs ---------------------------------------\n");
my (%pep_com, %pep_resc, %pro_resc);
%pep_com=%{clone(\%pep_acp)};
rescuePSMs($parahash{batches}, \%pep_acp, \%pep_fine, \%pep_com, \%pep_resc, \%pep2pro, \%prohash, \%pro_resc, \%scoreCutoff,$parahash{enable_group_specific_Jscore});

# [filter candidates] consolidate rescued PSMs: ensure one outfile maximally conveys only one peptide
consolidateRescuedPSMs(\%out2pep,\%pep_resc,\%pep_com, \%pro_resc, \%prohash,\%pep2pro);

# estimate FDR for rescued PSMs/peptides
#printLog("\nCalculate FDR for rescued peptides:\n");
#printFDR(\%pep_resc, \%pro_resc, $parahash{batches},0,1,0,0);

# estimate FDR for final (accepted + rescued) PSMs/peptides
printLog("\n--------------------------------------- After PSM rescue ---------------------------------------\n");
#printLog("\nCalculate FDR after merging originally accepted and rescued peptides:\n");
printFDR(\%pep_com, \%prohash, $parahash{batches},1,1,1,1,$summaryTable{batch_improved});

# part 4: output
printLog("\nOutputing results (by re-running -f internally; please be patient)\n");

printHash(\%pep_acp,"$parahash{output_folder}/intermediate/accepted_peptides.txt",$parahash{batches});
printHash(\%pep_resc,"$parahash{output_folder}/intermediate/rescued_peptides.txt",$parahash{batches});
printHash(\%pep_com,"$parahash{output_folder}/intermediate/final_combined_peptides.txt",$parahash{batches});

# stategy: foreach batch, mimic -s output based on %pep_com, and feed it to -f; -f will generate output folder for each batch, which can be used for downstream analysis (-l and -q)
foreach my $b (sort keys %{$parahash{batches}}) {
	# step 1: make the folder structure
	if (!(-e "$parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1")) {
		system(qq(mkdir -p $parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1));
	}

	my $sFolderPath=$parahash{paths}{pepXML};
	if ($searchEngine eq 'comet') {
		$sFolderPath = dirname($sFolderPath);
		system(qq(cp $sFolderPath/comet.params $parahash{output_folder}/$b/acceptedPSMs));
	} else { # jump or sequest
		$sFolderPath =~ s/\.pepXML$//;
		system(qq(cp $sFolderPath/jump.params $parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1));
	}
	#system(qq(ln -s ));

	# step 2: build %acceptedOutHash and printPepXML
	my (%acceptedOutHash);
	buildAcceptedPSMpepXML($b, \%pepXML_outhash, \%pep_com, \%pepXML_parahash, \%acceptedOutHash);

	my $pepxps=pepXML_parser->new;
	$pepxps->printPepXML($pepXML_parahash{$b},\%acceptedOutHash,"$parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1",1);

	# step 3: prepare -f .params
	# generate -f .params
	system(qq(cd $parahash{output_folder}/$b; jump -params));
	# configure -f .params
	my $dir = getcwd;

	if ($searchEngine eq 'comet') {
		params_f_update("$dir/$parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1","$parahash{output_folder}/$b/ParameterFiles/TMThh/jump_fc_HH_tmt10_human.params",'out',"$parahash{output_folder}/$b/jump_f.params",$parahash{mods},$parahash{pit_file},$parahash{database});
	} else { # jump or sequest
		params_f_update("$dir/$parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1","$parahash{output_folder}/$b/ParameterFiles/TMThh/jump_fj_HH_tmt10_human.params",'out',"$parahash{output_folder}/$b/jump_f.params",$parahash{mods},$parahash{pit_file},$parahash{database});
	}

	# step 4: run -f
	if (defined($parahash{jump_f_path}) and $parahash{jump_f_path} ne '0') {
		system(qq(cd $parahash{output_folder}/$b; perl $parahash{jump_f_path} jump_f.params));
	} else {
		# use program in production
		if (defined($parahash{dispatch}) and $parahash{dispatch} ne '0') {
			system(qq(cd $parahash{output_folder}/$b; jump -f jump_f.params --dispatch=$parahash{dispatch}));
		} else {
			system(qq(cd $parahash{output_folder}/$b; jump -f jump_f.params));
		}
	}

	# step 5: update ID.txt outfile path
	if ($parahash{mods} eq '0') {
		updateIDtxtPaths("$parahash{output_folder}/$b/sum_out/ID.txt",$pepxmlPath{$b},$searchEngine);
	} else {
		updateIDtxtPaths("$parahash{output_folder}/$b/sum_out_mod/IDmod.txt",$pepxmlPath{$b},$searchEngine);
	}
}

# print summary
printLog("\n------ Result Summary ------\n");
printLog("\nPeptide ID table:\n");
printSummaryTable(\%summaryTable,'peptide');
printLog("\n\nProtein ID table:\n");
printSummaryTable(\%summaryTable,'protein');

# finish
#system(qq(cat $parahash{output_folder}/jump_f_batch.log)); # re-fresh the summary
printLog("\njump_f_batch.pl is finished.\n\n");

#-------------------------------------------------------------------------------------
sub printSummaryTable {
# goal: print summary table with pretty format, hopefully
	my ($summaryTable,$level)=@_;
	my ($n1,$n2,$fdr1,$fdr2,@array1,@array2);

	printLog("\n\t\toriginal ID (FDR %)\tbatch-improved ID (FDR %)\n");
	printLog("----------------------------------------------------------\n");

	# print batch specific infor
	foreach my $batch (sort keys %{$$summaryTable{original}{$level}{batches}}) {

		$n1=$$summaryTable{original}{$level}{batches}{$batch}{n};
		$fdr1=$$summaryTable{original}{$level}{batches}{$batch}{FDR};
		$n2=$$summaryTable{batch_improved}{$level}{batches}{$batch}{n};
		$fdr2=$$summaryTable{batch_improved}{$level}{batches}{$batch}{FDR};

		printLog("$batch\t\t\t$n1 ($fdr1)\t\t$n2 ($fdr2)\n");

		push @array1,$n1;
		push @array2,$n2;
	}

	# print summary infor
	printLog("----------------------------------------------------------\n");

	$n1=$$summaryTable{original}{$level}{combine}{n};
	$fdr1=$$summaryTable{original}{$level}{combine}{FDR};
	$n2=$$summaryTable{batch_improved}{$level}{combine}{n};
	$fdr2=$$summaryTable{batch_improved}{$level}{combine}{FDR};
	printLog("combined\t\t$n1 ($fdr1)\t\t$n2 ($fdr2)\n");

	$n1=$$summaryTable{original}{$level}{overlap}{n};
	$fdr1=$$summaryTable{original}{$level}{overlap}{FDR};
	$n2=$$summaryTable{batch_improved}{$level}{overlap}{n};
	$fdr2=$$summaryTable{batch_improved}{$level}{overlap}{FDR};
	printLog("overlapped\t\t$n1 ($fdr1)\t\t$n2 ($fdr2)\n");

	#printLog("overlapped \%\t$p1\%\t\t\t$p2\%\n");
	printLog("overlapped \%:\n");

	my $p1=sprintf("%.2f",$n1*100/$$summaryTable{original}{$level}{combine}{n});
	my $p2=sprintf("%.2f",$n2*100/$$summaryTable{original}{$level}{combine}{n});
	printLog("  /(combined): \t\t$p1\%\t\t\t$p2\%\n");
	#printLog("  vs. combined: \t$p1\%\t\t\t$p2\%\n");

	$p1=sprintf("%.2f",$n1*100/mean(@array1));
	$p2=sprintf("%.2f",$n2*100/mean(@array1));
	#printLog("  vs. orig. batch mean: $p1\%\t\t\t$p2\%\n");
	printLog("  /(orig. batch mean):  $p1\%\t\t\t$p2\%\n");
}


sub updateIDtxtPaths {
# goal: update ID.txt outfile path
# procedures:
# 1) re-organize %pepxmlPath: %paths{$fraction}=$path;
# 2) parse ID.txt: loop for each row
#   a) parse $outfile and get $fraction
#   b) if exist %paths{$fraction}: update absolute path according to %pepxmlPath and print
	my ($IDtxt, $pepxmlPath, $searchEngine)=@_;

	# 0) backup file
	my $IDtxtbk="$IDtxt\_backup";
	system(qq(cp $IDtxt $IDtxtbk));

	# 1) re-organize %pepxmlPath: %paths{$fraction}=$path;
	my %paths;
	foreach my $p (keys %{$pepxmlPath}) {
		# /home/bbai/ADproject/human/individual/whole/batch2/w28/w28.1.pepXML
		my $file=basename($p);
		my $fraction=(split /\./,$file)[0];
		if ($searchEngine eq 'comet') {
			$p =~ s/\.pep.xml$//;
		} elsif ($searchEngine eq 'jump') { 
			$p =~ s/\.pepXML$//;
		} elsif ($searchEngine eq 'sequest') { 
			$p =~ s/\.pepcXML$//;
		}
		$paths{$fraction}=$p;
	}

	# 2) parse ID.txt: loop for each row
	open(IN,$IDtxtbk) or die "cannot open $IDtxtbk\n";
	open(OUT,">$IDtxt");

	while (<IN>) {
		chomp;
		if (/^Database/ or /^Peptide/) {
			print OUT "$_\n";
			next;
		}

		my @t=split /\;/;
		my $outfilepath=$t[2];

		my $outfile=basename($outfilepath);
		my $fraction=(split /\./,$outfile)[0];

		if (defined($paths{$fraction})) {
			$t[2]="$paths{$fraction}/$outfile";
			print OUT join(";",@t), "\n";
		} else {
			print "WARNING: no matching fraction ($fraction) for outfile: $outfilepath\n";
		}
	}

	close IN;
	close OUT;
}

sub params_f_update {
# goal: configure -f .params
# output folder
# -s path
# initial_outfile_fdr = 100
# FDR = 100
# bypass_filtering = 1
# mass_accuracy  = 0 # for mass calibration
# multistep_FDR_filtering = 0
# one_hit_wonders_removal = 0
# mods = $mods
# output_pepXML = 0
#
# usage example:
# params_f_update("$dir/$parahash{output_folder}/$b/acceptedPSMs/acceptedPSMs.1","$parahash{output_folder}/$b/ParameterFiles/TMThh/jump_fj_HH_tmt10_human.params",'out',"$parahash{output_folder}/$b/jump_f.params",$parahash{mods});

	my ($inputPath,$param,$outname,$output,$mods,$pitFile,$database)=@_;

	# read -f params template
	open(IN,$param) or die "Cannot open $param!!!\n";
        my @lines=<IN>;
        close IN;

	open(OUT,'>',$output);
        print OUT "$outname\: $inputPath\n"; # set input -s path

        for (my $i=0; $i<=$#lines;$i++) {
                $lines[$i] =~ s/^\s+//; # rm space at front
                next if ($lines[$i] =~ /^#/); # skip comment lines

                if ( $lines[$i]  =~ / = /  ) {
			my ($key, $value)=split / = /,$lines[$i];

			if ($key eq 'initial_outfile_fdr') {
				print OUT "$key = 100\n";
			} elsif ($key eq 'FDR') {
				print OUT "$key = 100\n";
			} elsif ($key eq 'one_hit_wonders_removal') {
				print OUT "$key = 0\n";
			} elsif ($key eq 'mods') {
				print OUT "$key = $mods\n";
			} elsif ($key eq 'pit_file') {
				print OUT "$key = $pitFile\n";
				print OUT "database = $database\n"; # temporary solution; should be revised if database is formally added to jump -params
			} elsif ($key eq 'output_pepXML') {
				print OUT "$key = 0\n";
			} elsif ($key eq 'mass_accuracy') {
				print OUT "$key = 0\n";
			} elsif ($key eq 'bypass_filtering') {
				print OUT "$key = 1\n";
			} elsif ($key eq 'multistep_FDR_filtering') {
				print OUT "$key = 0\n";
			} else {
				print OUT $lines[$i];
			}
		}
        }
	#print OUT "XCorr = 0\n"; # for mass calibration

        close OUT;

}

sub buildAcceptedPSMpepXML {
# goal: extract $outfile from %pepXML_outhash to %acceptedOutHash if $outfile is within %pep_com
	my ($b,$pepXML_outhash,$pep_com,$pepXML_parahash,$acceptedOutHash)=@_;

	my $pepxps=pepXML_parser->new;
	foreach my $pep (keys %{$pep_com}) {
		foreach my $outfile (keys %{$$pep_com{$pep}{$b}}) {
			my $hit=$$pep_com{$pep}{$b}{$outfile}{hit};

			$outfile =~ s/\.out$//;
			$outfile =~ s/\.spout$//;

			$pepxps->add_outfile_specificHit($$pepXML_parahash{$b}, $$pepXML_outhash{$b}, $acceptedOutHash, $outfile, $hit);
		}
	}
}

sub buildOutHash {
# goal: build %pepXML_outhash: add outfiles that pass Jscore cutoff
	my ($pepXML_outhash, $paraHash, $outInforHash, $min_Jscore, $jumpf_outfiles, $searchEngine)=@_;

	my $pepxps=pepXML_parser->new;
	foreach my $outfile (keys %{$outInforHash}) {
		next unless ($$outInforHash{$outfile}{hitShown}>0); # skip blank outfiles
		my $Jscore;
		if ($searchEngine eq 'jump') {
			$Jscore=$$outInforHash{$outfile}{WeightedEvalue}[1];
		} else { # comet or sequest
			$Jscore=$$outInforHash{$outfile}{xcorr}[1];
		}
		#next unless ($Jscore>=$min_Jscore); # skip low score outfiles
		#next unless ($Jscore>=$min_Jscore or scalar(keys %{$scanhash->{$outfile}})>0); # skip low score outfiles and those not already accepted
		next unless ($Jscore>=$min_Jscore or defined($$jumpf_outfiles{$outfile}) and $$jumpf_outfiles{$outfile}>0); # skip low score outfiles and those not already accepted

		$pepxps->add_outfile($paraHash, $outInforHash, $pepXML_outhash, $outfile);
	}
}

sub consolidateRescuedPSMs {
# goal: consolidate rescued PSMs: ensure one outfile maximally conveys only one peptide
# if multiple peptides are selected: select the top hit only; removed other hits from %pep_resc and %pep_com
# based on %out2pep
# loop for each outfile
# for: 1st hit to the end: if selected => remove all below records in %pep_resc and %pep_com
	my ($out2pep,$pep_resc,$pep_com,$pro_resc,$prohash,$pep2pro)=@_;

	foreach my $b (keys %{$out2pep}) {
		foreach my $outfile (keys %{$$out2pep{$b}}) {
			my $selected=0; # outfile selected in the rescue process?
			for (my $i=1; $i<scalar(@{$$out2pep{$b}{$outfile}}); $i++) {
				my $pep=$$out2pep{$b}{$outfile}[$i];
				#if (defined($$pep_resc{$pep}{$b}{$outfile})) {
				if (exists($$pep_resc{$pep}) 
				and exists($$pep_resc{$pep}{$b})
				and exists($$pep_resc{$pep}{$b}{$outfile})) {
					if ($selected==0) { # toppest hit to be selected
						$selected=1;
					} else { # one of the upper hit has been selected
						# remove this record in %pep_resc and %pep_com
						my $pro=$$pep2pro{$pep};

						# for rescued hashe
						my $batchRm=removeOutfileFromHash($pep_resc,$pep,$b,$outfile);

						if ($batchRm) {
							removeOutfileFromHash($pro_resc,$pro,$b,$pep);
						}

						# for combined hases
						$batchRm=removeOutfileFromHash($pep_com,$pep,$b,$outfile);
						if ($batchRm) {
							removeOutfileFromHash($prohash,$pro,$b,$pep);
						}
					}
				}
			}
		}
	}
}

sub removeOutfileFromHash {
# goal: recursively remove the most basic element (e.g., an outfile, a peptide) from hash
	my ($pephash,$pep,$b,$outfile)=@_;

	my $batchRm=0; # the outfiles for the whole batch are gone?
	delete $$pephash{$pep}{$b}{$outfile};
	if (scalar(keys %{$$pephash{$pep}{$b}})==0) {
		delete $$pephash{$pep}{$b};
		$batchRm=1;
		if (scalar(keys %{$$pephash{$pep}})==0) {
			delete $$pephash{$pep};
		}
	}

	return $batchRm;
}

sub printGroupScoreCuts {
# goal: print to temparory file: %groupScore{$group}
	my ($scoreCutoff,$output)=@_;

	open(OUT,">$output");
	print OUT "z\tpepLength\tscoreCut\n";
	foreach my $z (sort keys %{$scoreCutoff}) {
		foreach my $pepGroup (sort keys %{$$scoreCutoff{$z}}) {
			print OUT "$z\t$pepGroup\t$$scoreCutoff{$z}{$pepGroup}\n";
		}
	}
	close OUT;
}

sub groupScore {
# goal: build %scoreCutoff{$z}{$pepGroup}=$scoreCut
# based on %pep_acp
# $group= peptide2Group($pep)
# build push @{%scores{$group}}, $Jscore
# after looping all outfiles: %groupScore{$group}=quantile(@{%scores{$group}}, 0.95)
# print to temparory file: %groupScore{$group}
	my ($pep_acp,$scoreCutoff,$lowestScorePct)=@_;

	# group PSMs by charge states and peptide length, and recored scores in %scores
	my %scores;
	foreach my $pep (keys %{$pep_acp}) {

		my $pepGroup=pep2group($pep);
		foreach my $b (keys %{$$pep_acp{$pep}}) {
			foreach my $outfile (keys %{$$pep_acp{$pep}{$b}}) {
				# w28.34408.1.2.spout
				my $z=(split /\./,$outfile)[3];
				push @{$scores{$z}{$pepGroup}}, $$pep_acp{$pep}{$b}{$outfile}{Jscore};
			}
		}
	}

	# obtain lowest 5% (may set as parameter) quantile score for each group
	foreach my $z (keys %scores) {
		foreach my $pepGroup (keys %{$scores{$z}}) {
			$$scoreCutoff{$z}{$pepGroup}=quantile($lowestScorePct,@{$scores{$z}{$pepGroup}});
		}
	}
}

sub quantile {
# goal: within an array, return the value that correspond to a specific quantifle (0-1)
# currently, assume the quantile value corresponds to the ascendence of the array, i.e., 0% is the min, while 100% is max
	my ($q,@array)=@_;

	return '' unless (0<=$q and $q<=1);

	my $topPct=int(@array*$q);
	@array=sort {$a<=>$b} @array;

	return $array[$topPct];
}

sub pep2group {
# goals: assign a peptide to a group (defined by charge state and peptide length)
# charge states: 1-5 (not applicable here)
# peptide length: 7-20 (individual groups), 21-30, >30
# return: 7, 21-30., etc.
	my ($pep)=@_;

	my @range;
	$range[0]=20;
	$range[1]=30;

	my $l=length($pep);
	my $pepl='';

	if ($l<=$range[0]) {
		$pepl=$l;
	} elsif ($range[0]<$l and $l<=$range[1]) {
		$pepl=$range[0]+1;
		$pepl.="-$range[1]";
	} else {
		$pepl=$range[1]+1;
		$pepl.="+";
	}

	return $pepl;
}

sub proInfor {
# goals: parse id_uni_pep to establish relationship between peptides and proteins
# 1. build %pep2pro{$pep}=$protein
# 2. build %prohash{$pro}{$batch}{$pep}{TD}
	my ($id_uni_pep,$pep2pro,$prohash,$b)=@_;

	my %fpephash=%{retrieve($id_uni_pep)};

	foreach my $pep (keys %fpephash) {
		my $pro=$fpephash{$pep}{best_protein};
		my $TD=($pro =~ m/Decoy/)?'decoy':'target';

		$$pep2pro{$pep}=$pro;
		$$prohash{$pro}{$b}{$pep}{TD}=$TD;
	}
}

sub printLog {
# goal: print simualtaneously on screen and log file
	my ($sentence)=@_;

	open(LOG,">>$log_file");

	print $sentence;
	print LOG $sentence;

	close LOG;
}

sub printHash {
# goal: print hash of %pep_acp format
# output format: peptide, batch1: PSM#, outfiles, batch2: PSM#, outfiles etc.
	my ($pep_hash, $output, $batches)=@_;

	open(OUT,">$output");

	print OUT "peptide";
	foreach my $b (keys %{$batches}) {
		print OUT "\t$b.PSMs\t$b.outfiles";
	}
	print OUT "\tTD\n";
	
	foreach my $pep (keys %{$pep_hash}) {
		print OUT "$pep";
		my $TD='';
		foreach my $b (keys %{$batches}) {
			if (defined($$pep_hash{$pep}{$b})) {
				print OUT "\t",scalar(keys %{$$pep_hash{$pep}{$b}}),"\t";
				foreach my $outfile (keys %{$$pep_hash{$pep}{$b}}) {
					print OUT "$outfile,";
					$TD=$$pep_hash{$pep}{$b}{$outfile}{TD};
				}
			} else {
				print OUT "\t0\t";
			}
		}
		print OUT "\t$TD\n";
	}

	close OUT;
}

sub printFDR {
# goal: a wrapper to print FDR for 4 cases
	my ($pep_hash, $prohash, $batches, $printIndividual, $printUnion, $printIntersect, $printBS,$summaryTable)=@_;

	my $nBatch=scalar(keys %{$batches});

	# 4 cases:
	# print Individual batch
	if ($printIndividual) {
		printLog("\nPeptides that are identified in individual batch:\n");
		foreach my $b (keys %{$batches}) {
			printLog("\nPeptides that are identified in $b:\n");

			my ($target, $decoy, $FDR)=FDR_batch($pep_hash,$b);
			$$summaryTable{peptide}{batches}{$b}{n}=$target;
			$$summaryTable{peptide}{batches}{$b}{FDR}=$FDR;

			($target, $decoy, $FDR)=FDR_batch($prohash,$b,1);
			printLog "  Proteins: $target targets and $decoy decoys (FDR = $FDR%)\n";
			$$summaryTable{protein}{batches}{$b}{n}=$target;
			$$summaryTable{protein}{batches}{$b}{FDR}=$FDR;
		}

		printLog("\nPeptides that are identified after merging $nBatch batches:\n");
	}
	
	# overall (union)
	if ($printUnion) {
		printLog("\nPeptides that are identified in any batch (union):\n");

		my ($target, $decoy, $FDR)=FDR($pep_hash, 1, $nBatch);
		$$summaryTable{peptide}{combine}{n}=$target;
		$$summaryTable{peptide}{combine}{FDR}=$FDR;

		($target, $decoy, $FDR)=FDR($prohash, 1, $nBatch,1);
		printLog "  Proteins: $target targets and $decoy decoys (FDR = $FDR%)\n";
		$$summaryTable{protein}{combine}{n}=$target;
		$$summaryTable{protein}{combine}{FDR}=$FDR;
	}

	# overlapped (intersection)
	if ($printIntersect) {
		printLog("\nPeptides that are identified in all batches (intersection):\n");

		my ($target, $decoy, $FDR)=FDR($pep_hash, $nBatch, $nBatch);
		$$summaryTable{peptide}{overlap}{n}=$target;
		$$summaryTable{peptide}{overlap}{FDR}=$FDR;

		($target, $decoy, $FDR)=FDR($prohash, $nBatch, $nBatch,1);
		printLog "  Proteins: $target targets and $decoy decoys (FDR = $FDR%)\n";
		$$summaryTable{protein}{overlap}{n}=$target;
		$$summaryTable{protein}{overlap}{FDR}=$FDR;
	}

	# batch specific (union - intersect)
	if ($printBS) {
		printLog("\nPeptides that are identified in some batch(es) but not all (union - intersect):\n");
		FDR($pep_hash, 1, $nBatch-1);
		my ($target, $decoy, $FDR)=FDR($prohash, 1, $nBatch-1,1);
		printLog "  Proteins: $target targets and $decoy decoys (FDR = $FDR%)\n";
	}
}

sub FDR_batch {
# goals: calculate FDR for %pep_acp (only considering peptides that belong to a specific batch (as long as the batch has that peptide))

	my ($pep_acp, $batch, $silent)=@_;
	$silent ||= 0;

	my ($PSM_target, $PSM_decoy, $pep_target, $pep_decoy);
	$PSM_target=$PSM_decoy=$pep_target=$pep_decoy=0;
	foreach my $pep (keys %{$pep_acp}) {
		my $nBatch=scalar(keys %{$$pep_acp{$pep}});
		next unless (defined($$pep_acp{$pep}{$batch}));

		# count targets and decoys
		my $TD='';
		# PSM level
		foreach my $b (keys %{$$pep_acp{$pep}}) {
			next unless ($b eq $batch);
			foreach my $outfile (keys %{$$pep_acp{$pep}{$b}}) {
				$TD=$$pep_acp{$pep}{$b}{$outfile}{TD};
				if ($$pep_acp{$pep}{$b}{$outfile}{TD} eq 'target') {
					$PSM_target++;
				} else {
					$PSM_decoy++;
				}
			}
		}
		# peptide level
		if ($TD eq 'target') {
			$pep_target++;
		} else {
			$pep_decoy++;
		}
	}

	my $psmFDR=sprintf "%.2f",$PSM_decoy*100/($PSM_target+0.01);
	my $pepFDR=sprintf "%.2f",$pep_decoy*100/($pep_target+0.01);

	unless ($silent) {
		my $sentence=sprintf "  Outfiles: $PSM_target targets and $PSM_decoy decoys (FDR = %.2f%%)\n", $psmFDR;
		printLog($sentence);

		$sentence=sprintf "  Peptides: $pep_target targets and $pep_decoy decoys (FDR = %.2f%%)\n", $pepFDR;
		printLog($sentence);
	}

	return ($pep_target, $pep_decoy, $pepFDR);
}

sub rescuePSMs {
# goals: directly copy PSMs from %pep_fine (e.g., Jscore>20) to %pep_com (final peptide hash)
	my ($batches, $pep_acp, $pep_fine, $pep_com, $pep_resc, $pep2pro, $prohash, $pro_resc,$scoreCutoff,$enable_group_specific_Jscore)=@_;

	foreach my $pep (keys %{$pep_acp}) {
		# skip if all batch are fulfilled
		next if (scalar(keys %{$$pep_acp{$pep}})==scalar(keys %{$batches}));

		# if some batches are missing, loop for all the missing batches:
		foreach my $b (keys %{$batches}) {
		#foreach my $b (keys %{$pep_acp{$pep}}) {
			next if (defined($$pep_acp{$pep}{$b}));

			# if exists %pep_fine{$pep}{batch2}: add the highest score PSM into %pep_res{$pep}{batch2} and %pep_com{$pep}{batch2}
			# this function can be further revised to consider the neighbor peptide overlapped
			if (defined($$pep_fine{$pep}{$b})) {
				my $outfile=bestOutfile($$pep_fine{$pep}{$b},'Jscore');

				my $pepGrp=pep2group($pep);
				my $z=(split /\./,$outfile)[3];
				my $Jscore=$$pep_fine{$pep}{$b}{$outfile}{Jscore};

				# check if score > group specific cutoff
				if ($enable_group_specific_Jscore==0 or 
				defined($$scoreCutoff{$z}{$pepGrp}) and $Jscore>=$$scoreCutoff{$z}{$pepGrp}) { 
					$$pep_resc{$pep}{$b}{$outfile}=clone($$pep_fine{$pep}{$b}{$outfile});
					$$pep_com{$pep}{$b}{$outfile}=clone($$pep_fine{$pep}{$b}{$outfile});

					my $pro=$$pep2pro{$pep};
					my $TD=($pro =~ m/Decoy/)?'decoy':'target';
					$$prohash{$pro}{$b}{$pep}{TD}=$TD;
					$$pro_resc{$pro}{$b}{$pep}{TD}=$TD;
				}
			}

		}
	}
}

sub bestOutfile {
# Goal: pick up the best outfile
# criteria: select the one with max (min) of Jscore (or other field)
	#my ($hash, $field, $preferLarge)=@_;
	my ($hash, $field)=@_;

	my ($bestOut,$bestScore)=('',0);
	foreach my $outfile (keys %{$hash}) {
		if ($bestScore<$$hash{$outfile}{$field}) {
			$bestScore=$$hash{$outfile}{$field};
			$bestOut=$outfile;
		}
	}

	return $bestOut;
}

sub FDR {
# goals: calculate FDR for %pep_acp (only considering peptides that meet the range of number of available batches)

	my ($pep_acp, $minBatchN, $maxBatchN, $silent)=@_;
	$silent ||= 0;

	my ($PSM_target, $PSM_decoy, $pep_target, $pep_decoy);
	$PSM_target=$PSM_decoy=$pep_target=$pep_decoy=0;
	foreach my $pep (keys %{$pep_acp}) {
		# only considering peptides that meet the range of number of available batches
		my $nBatch=scalar(keys %{$$pep_acp{$pep}});
		next unless ($minBatchN<=$nBatch and $nBatch<=$maxBatchN);

		# count targets and decoys
		my $TD='';
		# PSM level
		foreach my $b (keys %{$$pep_acp{$pep}}) {
			foreach my $outfile (keys %{$$pep_acp{$pep}{$b}}) {
				$TD=$$pep_acp{$pep}{$b}{$outfile}{TD};
				if ($$pep_acp{$pep}{$b}{$outfile}{TD} eq 'target') {
					$PSM_target++;
				} else {
					$PSM_decoy++;
				}
			}
		}
		# peptide level
		if ($TD eq 'target') {
			$pep_target++;
		} else {
			$pep_decoy++;
		}
	}

	my $psmFDR=sprintf "%.2f",$PSM_decoy*100/($PSM_target+0.01);
	my $pepFDR=sprintf "%.2f",$pep_decoy*100/($pep_target+0.01);

	unless ($silent) {
		my $sentence=sprintf "  Outfiles: $PSM_target targets and $PSM_decoy decoys (FDR = %.2f%%)\n", $psmFDR;
		printLog($sentence);

		$sentence=sprintf "  Peptides: $pep_target targets and $pep_decoy decoys (FDR = %.2f%%)\n", $pepFDR;
		printLog($sentence);
	}

	return ($pep_target, $pep_decoy, $pepFDR);
}

sub pepxmlOuthash2pephash {
# goals: 
# 1. add this batch specific hash %outInforHash to %pep_fine
# 2. convert outfile hash to peptide hash: %pep_fine{$pep}{batch1}{$outfile}{Jscore/nbPep/TD/scan/fraction}

	my ($batch, $paraHash, $outInforHash, $pep_fine, $scoreCutoff, $multiHit_max_dJn, $out2pep, $searchEngine)=@_;

	foreach my $outfile (keys %{$outInforHash}) {
		next unless ($$outInforHash{$outfile}{hitShown}>0); # skip blank outfiles

		# collect outfile and peptide information, with two requirements:
		# 1. Jscore > threshold
		# 2. if close scores for 2nd and below hits (determined by dJn), these hits are also considered
		for (my $i=1; $i<=$$outInforHash{$outfile}{hitShown}; $i++) {
			my $dJn=$$outInforHash{$outfile}{deltacn}[$i];
			
			#my $Jscore=$$outInforHash{$outfile}{WeightedEvalue}[$i];
			my $Jscore;
			if ($searchEngine eq 'jump') {
				$Jscore=$$outInforHash{$outfile}{WeightedEvalue}[1];
			} else { # comet or sequest
				$Jscore=$$outInforHash{$outfile}{xcorr}[1];
			}
			last unless ($Jscore>=$scoreCutoff); # skip low score PSMs

			my ($fraction, $scan, $ppi, $z)=split /\./,$outfile;
			#my $pep=$$outInforHash{$outfile}{peptide}[$i];
			my $pep=$$outInforHash{$outfile}{fullPeptide}[$i];
			my $pro=$$outInforHash{$outfile}{protein}[$i];

			my $TD=($pro =~ m/Decoy/)?'decoy':'target';

			# build %pep_fine
			$$pep_fine{$pep}{$batch}{$outfile}{Jscore}=$Jscore;
			$$pep_fine{$pep}{$batch}{$outfile}{TD}=$TD;
			$$pep_fine{$pep}{$batch}{$outfile}{scan}=$scan;
			$$pep_fine{$pep}{$batch}{$outfile}{fraction}=$fraction;
			$$pep_fine{$pep}{$batch}{$outfile}{hit}=$i;

			# build %out2pep
			$$out2pep{$batch}{$outfile}[$i]=$pep;

			last if ($dJn>$multiHit_max_dJn);
		}
	}
}

sub add_idtxtPephash2pephash {
# goals: 
# 1. add batch specific pephash to %pep_acp
# 2. re-format: %pep_acp{$pep}{batch1}{$outfile}{Jscore/nbPep/TD/scan/fraction}

	my ($batch, $bs_pep_hash, $bs_scan_hash, $pep_acp)=@_;

	foreach my $pep (keys %$bs_pep_hash) {

		my $pro=(keys %{$bs_pep_hash->{$pep}->{proteins}})[0];
		my $TD=($pro =~ m/Decoy/)?'decoy':'target';

		foreach my $outfile (keys %{$bs_pep_hash->{$pep}->{outfiles}}) {
			my $fullOut=$outfile;
			# p39.83007.1.3.spout
			$outfile=basename($outfile);
			#my ($fraction, $scan, $ppi, $z)=split /\./,basename($outfile);
			my ($fraction, $scan, $ppi, $z)=split /\./,$outfile;

			$$pep_acp{$pep}{$batch}{$outfile}{Jscore}=$bs_scan_hash->{$fullOut}->{XCorr};
			$$pep_acp{$pep}{$batch}{$outfile}{TD}=$TD;
			$$pep_acp{$pep}{$batch}{$outfile}{fraction}=$fraction;
			$$pep_acp{$pep}{$batch}{$outfile}{scan}=$scan;
			$$pep_acp{$pep}{$batch}{$outfile}{hit}=1;
		}
	}
}

# basic statistics functions cp from http://www.perlmonks.org/?node_id=1089988
sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}


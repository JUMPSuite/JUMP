#!/usr/bin/perl

package Utils::PublicationTable;

use strict;
use warnings;
use Cwd;
use Excel::Writer::XLSX;
use File::Basename;
use Spreadsheet::XLSX;
use Utils::Parse;

sub new {
	my ($class) = @_;
	my $self = {};
	bless ($self,$class);
	return $self;
}

sub generatePublicationTables {
	my $self = shift @_;
	my ($params, $jumpfTable, $scanZeroFile, $publicationDir) = @_;
	my $idTxt = $$params{'idtxt'};
	my $idJumpf;
	if ($jumpfTable =~ /site/) {
		my $intermediateDir = getcwd() . "\/$$params{'save_dir'}" . "/intermediate/";
		$idJumpf = $intermediateDir . $jumpfTable;
	} else {
		$idJumpf = dirname($idTxt) . "/publications/". $jumpfTable;  # JUMP -f publication table
	}
	my $parse = Utils::Parse -> new();
	my ($idHeader, $idHash, $idAccArray) = $parse -> parseJumpfPublicationTable($idJumpf);
	my ($normReporters, $normHash) = $parse -> parseNormalizedReporterIntensity($params, $jumpfTable);
	my $isComparison = $$params{'comparison_analysis'};
	my $compHash;
	if ($isComparison) {
		$compHash = $parse -> parseComparisonResult($params, $jumpfTable);
	}	

	## Common accessions from three hashes (i.e. idHash, normHash and compHash)
	my $commonAccs;
	map {$commonAccs -> {$_}++} keys %$idHash;
	map {$commonAccs -> {$_}++} keys %$normHash;
	map {$commonAccs -> {$_}++} keys %$compHash;
	my %commonAccs;
	if ($isComparison) {
		%commonAccs = map {($_, 1)} grep {$commonAccs -> {$_} > 2} keys %$commonAccs;
	} else {
		%commonAccs = map {($_, 1)} grep {$commonAccs -> {$_} > 1} keys %$commonAccs;
	}
	my @sortedAccs = map {$_ -> [1]} grep {exists $commonAccs{$_ -> [1]}} @$idAccArray;

	my $nQuantified = scalar(@sortedAccs);
	my $header = "";
	if ($jumpfTable eq "id_uni_prot.txt") {
		$header = "Unique proteins quantified by mass spectrometry (n = $nQuantified)";
	} elsif ($jumpfTable eq "id_all_prot.txt") {
		$header = "All proteins quantified by mass spectrometry (n = $nQuantified)";
	} elsif ($jumpfTable eq "id_uni_pep.txt") {
		$header = "Unique peptides quantified by mass spectrometry";
	} elsif ($jumpfTable eq "id_all_pep.txt") {
		$header = "All peptides quantified by mass spectrometry";
	} elsif ($jumpfTable eq "id_uni_site.txt") {
		$header = "Modification sites of unique proteins quantified by mass spectrometry (n = $nQuantified)";
	} elsif ($jumpfTable eq "id_all_site.txt") {
		$header = "Modification sites of all proteins quantified by mass spectrometry (n = $nQuantified)";
	} else {
		die "  jump -f publication table should be correctly specified\n";
	}
	
	## Peptide or protein publication tables?
	my $pepTable = 0;
	if ($jumpfTable =~ /pep/) {
		$pepTable = 1;
	}
	
	######################################################
	## Write tables										##
	## 1. id_uni/all_pep/prot/site_quan.xlsx			##
	## 2. id_uni/all_pep/prot/site_quan_with_na.xlsx	##
	## 3. id_uni/all_pep/prot/site_quan.txt				##
	######################################################
	## Open files
	## 1. PublicationTable
	my $table = $jumpfTable;
	$table =~ s/\.txt/_quan\.xlsx/;
	$table = $publicationDir."/".$table;
	my $workbook = Excel::Writer::XLSX -> new($table);
	my $tmpdir = File::Temp->newdir();
	$workbook -> set_tempdir($tmpdir->dirname);
	my $worksheet = $workbook -> add_worksheet();
	my $format = $workbook -> add_format();
	## 2. PublicationTableNA
	my $tableNA = $table;
	$tableNA =~ s/\.xlsx/_with_na\.xlsx/;
	my $workbookNA = Excel::Writer::XLSX -> new($tableNA);
	$workbookNA -> set_tempdir($tmpdir->dirname);
	my $worksheetNA = $workbookNA -> add_worksheet();
	my $formatNA = $workbookNA -> add_format();	
	## 3. PublicationTableTxt
	my $tableTxt = $table;
	$tableTxt =~ s/xlsx/txt/;
	
	## Write a header
	## 1. To publicationTable
	$format -> set_bold(1);
	$worksheet -> write_string(0, 0, $header, $format);
	## 2. To publicationTableNA
	$formatNA -> set_bold(1);
	$worksheetNA -> write_string(0, 0, $header, $formatNA);
	## 3. To publicationTableTxt
	open (TXT, ">", $tableTxt);
	print TXT "$header\n";
	## Additional header for the peptide publication tables
	if ($pepTable) {
		$header = "n = $nQuantified peptides";
		$worksheet -> write_string(1, 0, $header, $format);
		$worksheetNA -> write_string(1, 0, $header, $format);
	}
	
	## Write colunm names
	my @colNames = split(/\t/, $idHeader);
	push (@colNames, @{$normReporters});
	if ($isComparison) {
		foreach my $compInd (sort {$a <=> $b} keys %{$$params{'comparisonNames'}}) {
			my $comparisonName = $$params{'comparisonNames'}{$compInd};
			push (@colNames, @{$$compHash{'header'}{$comparisonName}});
		}
	}
	my $rowInd = 1;
	if ($pepTable) {
		$rowInd = 4;
	}
	map {$worksheet -> write_string($rowInd, $_, $colNames[$_], $format)} (0 .. $#colNames);	# publicationTable	
	map {$worksheetNA -> write_string($rowInd, $_, $colNames[$_], $formatNA)} (0 .. $#colNames);	# publicationTableNA

	## for jump -i, column names should only contain reporter names in .txt files
	my @normReportersTxt;
	for (my $i = 0; $i < scalar(@{$normReporters}); $i++) {
		my ($normReporter) = ($$normReporters[$i] =~ /(sig[0-9NC]+)\s/);
		push (@normReportersTxt, $normReporter);
	}
	print TXT "$idHeader\t", join("\t", @normReportersTxt), "\n";	# publicationTableTxt
	
	## Write values (i.e. information in ID(mod).txt, normalized intensities and statistical testing results)
	my $offset = $rowInd + 1;	
	my %used;
	my %protCountHash;	## for peptide publication tables
	my $protInd = 0;	## for publication tables of 'id_all_pep.txt'
	for (my $i = 0; $i < scalar(@sortedAccs); $i++) {	
		my $acc = $sortedAccs[$i];
		$used{$acc} = 1;
		my @values;
		## Peptide publication tables
		if ($pepTable && $jumpfTable eq "id_uni_pep.txt") {
			my $prot = (keys %{$$idHash{$acc}})[0];
			@values = (@{$$idHash{$acc}{$prot}}, @{$$normHash{$acc}}); 
			my $protGroup = (split(/\./, $$idHash{$acc}{$prot}[1]))[0];
			$protCountHash{$protGroup} = 1;
		} elsif ($pepTable && $jumpfTable eq "id_all_pep.txt") {
			if ($i > 0 && $acc eq $sortedAccs[$i - 1]) {
				$protInd++;
			} else {
				$protInd = 0;
			}
			my $prot = (sort {$$idHash{$acc}{$a}[1] cmp $$idHash{$acc}{$b}[1]} keys %{$$idHash{$acc}})[$protInd];
			@values = (@{$$idHash{$acc}{$prot}}, @{$$normHash{$acc}});
			$protCountHash{$prot} = 1;
		## Protein publication tables
		} else {
			@values = (@{$$idHash{$acc}}, @{$$normHash{$acc}});
		} 
		print TXT join("\t", @values), "\n";	# publicationTableTxt
		if ($isComparison) {
			foreach my $compInd (sort {$a <=> $b} keys %{$$params{'comparisonNames'}}) {
				my $comparisonName = $$params{'comparisonNames'}{$compInd};
				push (@values, @{$$compHash{$acc}{$comparisonName}});
			}			
		}
		map {$worksheet -> write($i + $offset, $_, $values[$_])} (0 .. $#values);	# publicationTable
		map {$worksheetNA -> write($i + $offset, $_, $values[$_])} (0 .. $#values);	# publicationTableNA	
	}
	close (TXT);
	
	## Additional header for peptide publication tables
	if ($pepTable) {
		my $nProts = scalar(keys %protCountHash);
		if ($jumpfTable eq "id_uni_pep.txt") {
			$header = "n = $nProts protein groups";
		} elsif ($jumpfTable eq "id_all_pep.txt") {
			$header = "n = $nProts proteins";
		} else {
			die "  jump -f table should be correctly specified; either id_uni_pep.txt or id_all_pep.txt\n";
		}		
		$worksheet -> write_string(2, 0, $header, $format);
		$worksheetNA -> write_string(2, 0, $header, $format);
	}
	$workbook -> close();
	
	## Write "NA"s to publicationTableNA
	my $i = 0;
	foreach my $acc (keys %$idHash) {
		next if ($used{$acc});
		my @values;
		if ($pepTable && $jumpfTable eq "id_all_pep.txt") {
			foreach my $prot (sort {$$idHash{$acc}{$a}[1] cmp $$idHash{$acc}{$b}[1]} keys %{$$idHash{$acc}}) {
				@values =  @{$$idHash{$acc}{$prot}};
				my @NAs = ("NA") x (scalar(@colNames) - scalar(@values));
				@values = (@values, @NAs);
				map {$worksheetNA -> write(scalar(@sortedAccs) + $i + $offset, $_, $values[$_])} (0 .. $#values);
				$i++;
			}
		} elsif ($pepTable && $jumpfTable eq "id_uni_pep.txt") {
			my $prot = (keys %{$$idHash{$acc}})[0];
			@values = @{$$idHash{$acc}{$prot}};
			my @NAs = ("NA") x (scalar(@colNames) - scalar(@values));
			@values = (@values, @NAs);
			map {$worksheetNA -> write(scalar(@sortedAccs) + $i + $offset, $_, $values[$_])} (0 .. $#values);
			$i++;
		} else {
			@values = @{$$idHash{$acc}};
			my @NAs = ("NA") x (scalar(@colNames) - scalar(@values));
			@values = (@values, @NAs);
			map {$worksheetNA -> write(scalar(@sortedAccs) + $i + $offset, $_, $values[$_])} (0 .. $#values);
			$i++;
		} 
	}
	$workbookNA -> close();
	
	######################################################
	## Write tables										##
	## 4. id_uni/all_pep/prot_zero_psm_quan.xlsx		##
	## 5. id_uni/all_pep/prot_zero_psm_only_quan.xlsx	##
	## No zero- tables for site-level					##
	######################################################
	if ($jumpfTable =~ /site/) {
		return;
	}
	## Open files
	## 4. PublicationTableFiltered
	my $tableFiltered = $table;
	$tableFiltered =~ s/quan\.xlsx/zero\_psm\_quan\.xlsx/;
	my $workbookFiltered = Excel::Writer::XLSX -> new($tableFiltered);
	$workbookFiltered -> set_tempdir($tmpdir->dirname);
	my $worksheetFiltered = $workbookFiltered -> add_worksheet();
	my $formatFiltered = $workbookFiltered -> add_format();
	## 5. PublicationTableFilteredOnly
	my $tableFilteredOnly = $table;
	$tableFilteredOnly =~ s/quan\.xlsx/zero\_psm\_only\_quan\.xlsx/;
	my $workbookFilteredOnly = Excel::Writer::XLSX -> new($tableFilteredOnly);
	$workbookFilteredOnly -> set_tempdir($tmpdir->dirname);
	my $worksheetFilteredOnly = $workbookFilteredOnly -> add_worksheet();
	my $formatFilteredOnly = $workbookFilteredOnly -> add_format();
	
	## Create a hash containing all the information in the "xxx_scan_zero.txt" file
	my %zeroInfo;
	open (ZERO, "<", $scanZeroFile) || die "  Cannot open $scanZeroFile\n";
	<ZERO>;
	<ZERO>;	## skip the first two lines
	while (<ZERO>) {
		chomp;
		my @elems = split(/;/, $_);
		## For the protein-level zero table
		next if (!defined $elems[1]);
		my ($pep, $prot, $outfile) = @elems[0, 1, 2];
		if ($pepTable) {
			$zeroInfo{$pep}{$prot}{$outfile} = $_;
		} else {
			$zeroInfo{$prot}{$pep}{$outfile} = $_;
		}
	}
	close (ZERO);
	
	## Write column names
	my $nReporters = scalar(@{$normReporters});
	my $nCompHeaders = 0;
	my @zeroColNames = split(/\t/, $idHeader);
	push (@zeroColNames, "PSM_scan");
	push (@zeroColNames, @{$normReporters});
	if ($isComparison) {
		foreach my $compInd (sort {$a <=> $b} keys %{$$params{'comparisonNames'}}) {
			my $comparisonName = $$params{'comparisonNames'}{$compInd};
			push (@zeroColNames, @{$$compHash{'header'}{$comparisonName}});
			$nCompHeaders += scalar(@{$$compHash{'header'}{$comparisonName}});
		}
	}
	unshift (@zeroColNames, "Zero_only");
	$formatFiltered -> set_bold(1);
	$formatFilteredOnly -> set_bold(1);
	map {$worksheetFiltered -> write_string(4, $_, $zeroColNames[$_], $formatFiltered)} (0..$#zeroColNames);
	map {$worksheetFilteredOnly -> write_string(4, $_, $zeroColNames[$_], $formatFilteredOnly)} (0..$#zeroColNames);
	
	## Write values
	$offset = 5;
	$i = 0;
	my $j = 0;
	my $zeroOnly = 0;	 
	undef %used;
	foreach my $idKey (keys %{$idHash}) {
		foreach my $zeroKey (keys %{$zeroInfo{$idKey}}) {
			foreach my $outfile (keys %{$zeroInfo{$idKey}{$zeroKey}}) {
				my @zeroInfoArray = split(/;/, $zeroInfo{$idKey}{$zeroKey}{$outfile});
				my @values;
				if ($pepTable) {
					my $prot = (keys %{$$idHash{$idKey}})[0];
					@values = (@{$$idHash{$idKey}{$prot}}, @zeroInfoArray[2, -$nReporters..-1]);
				} else {
					@values = (@{$$idHash{$idKey}}, @zeroInfoArray[2, -$nReporters..-1]);
				}
				my @NAs = ("NA") x $nCompHeaders;
				@values = (@values, @NAs);
				if (!defined $$compHash{$idKey}) {	## zero protein entry has not been quantified nor used for the comparison analysis
					$zeroOnly++;
					$worksheetFiltered -> write($i + $offset, 0, "zero_only");
					map {$worksheetFiltered -> write($i + $offset, 1 + $_, $values[$_])} (0 .. $#values);
					$i++;
					$worksheetFilteredOnly -> write($j + $offset, 0, "zero_only");
					map {$worksheetFilteredOnly -> write($j + $offset, 1 + $_, $values[$_])} (0 .. $#values);
					$j++;
				} else {	## zero protein entry has been quantified and used for the comparison study (by other peptides/PSMs)
					$worksheetFiltered -> write($i + $offset, 0, "mixed_zero");
					map {$worksheetFiltered -> write($i + $offset, 1 + $_, $values[$_])} (0 .. $#values);
					$i++;					
				}
				if (defined $$compHash{$idKey} && !defined $used{$idKey}) {
					my @values;
					if ($pepTable) {
						my $prot = (keys %{$$idHash{$idKey}})[0];
						@values = (@{$$idHash{$idKey}{$prot}}, "NA", @{$$normHash{$idKey}});
					} else {
						@values = (@{$$idHash{$idKey}}, "NA", @{$$normHash{$idKey}});
					} 
					foreach my $compInd (sort {$a <=> $b} keys %{$$params{'comparisonNames'}}) {
						my $comparisonName = $$params{'comparisonNames'}{$compInd};
						push (@values, @{$$compHash{$idKey}{$comparisonName}});
					}
					$worksheetFiltered -> write($i + $offset, 0, "mixed_non_zero");
					map {$worksheetFiltered -> write($i + $offset, 1 + $_, $values[$_])} (0 .. $#values);
					$i++;
					$used{$idKey} = 1;
				}
			}
		}
	}
	my $total = $i;
	my $prefix;
	if ($jumpfTable eq "id_uni_prot.txt") {
		$prefix = "Unique proteins ";
	} elsif ($jumpfTable eq "id_all_prot.txt") {
		$prefix = "All proteins ";
	} elsif ($jumpfTable eq "id_uni_pep.txt") {
		$prefix = "Unique peptides ";
	} elsif ($jumpfTable eq "id_all_pep.txt") {
		$prefix = "All peptides ";
	} else {
		die "  jump -f table should be correctly specified\n";
	}
	my $filteredHeader = $prefix."matched by TMT PSMs with zero intensity (n = $total)";
	$worksheetFiltered -> write_string(0, 0, $filteredHeader, $formatFiltered);
	$worksheetFilteredOnly -> write_string(0, 0, $filteredHeader, $formatFilteredOnly);
	my $nonZeroOnly = $total - $zeroOnly;
	$filteredHeader = $prefix."matched by both zero and non-zero intensity TMT PSMs (n = $nonZeroOnly)";
	$worksheetFiltered -> write_string(1, 0, $filteredHeader, $formatFiltered);
	$worksheetFilteredOnly -> write_string(1, 0, $filteredHeader, $formatFilteredOnly);
	$filteredHeader = $prefix."matched by only zero intensity TMT PSMs (n = $zeroOnly)";
	$worksheetFiltered -> write_string(2, 0, $filteredHeader, $formatFiltered);
	$worksheetFilteredOnly -> write_string(2, 0, $filteredHeader, $formatFilteredOnly);
	$workbookFiltered -> close();
	$workbookFilteredOnly -> close();
}

sub generateCombinedPublicationTable {
	my $self = shift @_;
	my ($params, $publicationDir, $protFile) = @_;
	my $pepFile = $publicationDir."/id_uni_pep_quan.xlsx";
	my $combinedFile = $publicationDir."/id_uni_pep_prot_combined_quan.xlsx";
	
	my @pepHeader;
	my $pepData;
	my $nPhosphopeptides = 0;
	my $excel = Spreadsheet::XLSX -> new ($pepFile);
	foreach my $sheet (@{$excel -> {Worksheet}}) {
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		foreach my $headerCol ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
			push(@pepHeader,$sheet -> {Cells} [4] [$headerCol]->{Val});
		}
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}) {
			$sheet -> {MaxCol} ||= $sheet -> {MinCol};
			my @row_data = ();
			foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
				my $cell = $sheet -> {Cells} [$row] [$col];
				if ($cell) {
					push (@row_data,$cell -> {Val});
				}
			}
			next if (!defined($sheet -> {Cells} [$row] [2]->{Val}));
			next if ($sheet -> {Cells} [$row] [0] eq "Peptides");
			$pepData->{$sheet -> {Cells} [$row] [0]->{Val}} = \@row_data;			
		}
		$nPhosphopeptides = $sheet -> {MaxRow}-1;
	}
	
	my @protHeader;
	my $protData;
	my $nPhosphoproteins = 0;
	$excel = Spreadsheet::XLSX -> new ($protFile);
	foreach my $sheet (@{$excel -> {Worksheet}}) {
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		foreach my $headerCol ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
			push(@protHeader,$sheet -> {Cells} [1] [$headerCol]->{Val});
		}
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}) {
			$sheet -> {MaxCol} ||= $sheet -> {MinCol};
			my @row_data = ();
			foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
				my $cell = $sheet -> {Cells} [$row] [$col];
				if ($cell) {
					push (@row_data,$cell -> {Val});
				}
			}
			next if(!defined($sheet -> {Cells} [$row] [1]->{Val}));
			$protData->{$sheet -> {Cells} [$row] [1]->{Val}} = \@row_data;
		}
		$nPhosphoproteins = $sheet -> {MaxRow}-1;
	}
	
	my $workbook = Excel::Writer::XLSX->new($combinedFile);
	my $tmpdir = File::Temp->newdir();
	$workbook->set_tempdir( $tmpdir->dirname );
	my $worksheet = $workbook->add_worksheet();
	my $format = $workbook->add_format();
	my $i = 4;
	my @colNames = (@pepHeader, @protHeader);
	$format->set_bold(1);
	map {$worksheet -> write_string(3, $_, $colNames[$_], $format)} (0 .. $#colNames);
	my $nPeptides = 0;
	my %protCount;		  
	foreach my $peptide (keys %$pepData) {
		my $protAcc = $$pepData{$peptide}[2];
		next if(!defined($protData->{$protAcc}));
		$nPeptides++;
		$protCount{$protAcc} = 1;
		my @pepArray = @{$pepData->{$peptide}};	
		my @protArray = @{$protData->{$protAcc}};
		my @values = (@pepArray, @protArray);
		map {$worksheet -> write( $i , $_, $values[$_])} (0 .. $#values);
		$i++;
	}
	my $nProteins = scalar(keys %protCount);
	my $header = "Unique peptides quantified by mass spectrometry";
	$format->set_bold(1);
	$worksheet->write_string(0, 0, $header, $format);
	$header = "n = $nPeptides peptides";
	$worksheet->write_string(1, 0, $header, $format);
	$header = "n = $nProteins proteins quantified by whole proteome analysis";
	$worksheet->write_string(2, 0, $header, $format);
	$workbook -> close();
}

1;

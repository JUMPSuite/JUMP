#!/usr/bin/perl 

##################################################
## Generation of JUMP/SEQUEST search databases	##
## Generation of PIT (protein inference table)	##
##################################################
use strict;
use warnings;
use Cwd;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin";
use Utils::DatabaseUtils;
use Utils::DatabaseGeneration;

######################
## Initialization	##
######################
if (scalar(@ARGV) != 1) {
	print "USAGE:\n\tjump -d jump_d.params\n";
	exit;
}
open (my $log, ">", "jump_d.log");
my $paramsFile = shift (@ARGV);
unless (-e $paramsFile) {
	print "Cannot open $paramsFile\n";
	print "You need to have builddb.params file in the current working directory\n";
	print $log "Cannot open $paramsFile\n";
	print $log "You need to have builddb.params file in the current working directory\n";
	exit;
}
print "jump -d initialized\n\n";
print $log "jump -d initialized\n\n";

##################################################
## Read builddb.params file and Initialize	##
##################################################
my $params = parseParams($paramsFile);
if (!-e $$params{'jump.params'}) {
	print "Cannot open $$params{'jump.params'}\n";
	print "Please check the path for \"jump.params\" in the builddb.params file carefully\n";
	print $log "Cannot open $$params{'jump.params'}\n";
	print $log "Please check the path for \"jump.params\" in the builddb.params file carefully\n";
	exit;
}
my $jumpParams = parseParams($$params{'jump.params'});
my $utils = Utils::DatabaseUtils -> new();
my $db = Utils::DatabaseGeneration -> new();

##################################################
## Automatic generation of a new database name	##
##################################################
my $dbSuffix = $utils -> generateDbName($jumpParams, $log);
my $dbName = $$params{'output_prefix'}."_".$dbSuffix;
print "\n";
print $log "\n";

##################################################################################
## Read input .fasta file(s) and 						##
## generate a new .fasta file and the corresponding .pit file for the database	##
## Contaminants and decoys will be added, if necessary				##
##################################################################################
$utils -> generateFastaAndPit($params, $dbName, $log);
print "\n";
print $log "\n";

##########################################################
## Generate a new JUMP or SEQUEST search database	##
##########################################################
if ($$params{'bypass_db_generation'} == 1) {
	print "Skipped the database generation (only .fasta and .pit files are generated)\n";
} else {
	$db -> generateDb($jumpParams, $dbName, $log);
	print "\n";
	print $log "\n";
}
close ($log);

##################
## Subroutine	##
##################
sub parseParams {
	my ($paramFile) = @_;
	my %params;
	open (FILE, "<", $paramFile) or die "Cannot open a parameter file\n";
	while (<FILE>) {
		if ($_ =~ s/\s*([;\#].*)$// ) {
			my $dummy = $1;
		}
		if ($_ =~ /^(.+?)\s*=\s*(.+)$/) {
			my ($key, $value) = ($1, $2);
			$value =~ s/\s+$//o;
			$params{$key} = $value;
		}
	}
	close (FILE);
	return (\%params);
}

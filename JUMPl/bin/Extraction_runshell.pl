use Getopt::Long;

use lib $ENV{"JUMP_L_LIB"};
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Spiders::ProcessingMzXML;
use Spiders::Params;
use Spiders::MassCorrection;

GetOptions('IDmod=s'=>\$IDmod,
			'outdir=s'=>\$outdir,
			'fraction=s'=>\$fraction,
			'parameter=s'=>\$parameter,
		);
		
	
my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();
my $frac_scan = read_IDmod($IDmod,$fraction);

print "  Extracting MS2 peaks from .mzXML\n";
#print LOGFILE "  Extracting MS2 peaks from .mzXML\n";
my $proc_xml = new Spiders::ProcessingMzXML();
my $basename = basename($fraction);
# my $new_path = getcwd() . "/$basename";
my $new_path = $outdir."/$basename";

$proc_xml ->set_dta_path($new_path);
my $mzXML = basename($fraction);
$mzXML =~s/(\.\d+)//;
$mzXML = dirname($fraction) . "/". $mzXML . ".mzXML";

$proc_xml ->set_mzXML_file($mzXML);

############### preprocessing #########################
my (%ms_hash,%msms_hash,@mz_array);

$proc_xml ->set_parameter($params);
$proc_xml->generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);	
my $masscorr = new Spiders::MassCorrection();
my ($msms_hash_corrected,$mz_array_corrected) = $masscorr->massCorrection(\%ms_hash, \%msms_hash, \@mz_array, $params);
%msms_hash = %$msms_hash_corrected;

foreach my $scan (keys %{$frac_scan->{$fraction}})
{
	my $mz = $frac_scan->{$fraction}->{$scan}->{'MH1'};
	my $scan_num = $scan;
	my @data = split(/\./,$scan_num);

	my $mz_array = $msms_hash{$data[1]}{'msms_mz'};
	my $int_array = $msms_hash{$data[1]}{'msms_int'};
	
	$proc_xml->generate_dta_file($scan,$mz,$mz_array,$int_array);
}

sub read_IDmod
{
	##########################################
	## Read ID.txt file and generate hashes ##
	##########################################
	my ($idTxt,$fraction) = @_;
	
	my %frac_scan;
	
	print "  Parsing information from $idTxt\n";
#	print LOGFILE "Parsing information from $idTxt\n";

	open (INFILE, $idTxt);
	while (<INFILE>) {
		chomp;
		## Header line
		## Finally, the header includes the header of ID.txt and ASCORE header (defined above) 
		if ($_ =~ /^Peptide/) {
			$_ =~ s/\;/\,/g;
			$header = $_.",".$header;
		## Database line
		} elsif ($_ =~ /^Database/) {
			$databaseHeader = $_;
		## Peptide information line
		} else {
			my @elems = split(/\;/, $_);
			if (($elems[0] ne "") && (defined($elems[0]))) {
				#####################
				## Generate hashes ##
				#####################
				my $line = $_;
				## Key: a directory name where .out files in ID.txt file come from (used to generate and indicate ACORE.csv files)
				## Value: .out file names comes from the "key" directory
				## For example, if ID.txt file contains .out files from  three fractions, 
				## %dir2out will have three keys (directory names) representing those fractions  
				my $dir = dirname($elems[2]);
				next if ($dir ne $fraction);
				my $file = basename($elems[2]);
				$file =~ s/(.out|.spout)//;

				
				$frac_scan{$dir}{$file}{'peptide'}=$elems[0];
				$frac_scan{$dir}{$file}{'MH1'}=$elems[4];
			
				## %out2id
				## Key: .out file name (with an absolute path)
				## Value: corresponding output information in id.txt file
				## Therefore, %out2id has the keys equal to the number of unique .out files in ID.txt file
				my $idKey = $dir."/".$file;
				$_ =~ s/\;/\,/g;
				push (@{$out2id{$idKey}}, $_);
			}
		}
	}
	close INFILE;
	return (\%frac_scan);
}	

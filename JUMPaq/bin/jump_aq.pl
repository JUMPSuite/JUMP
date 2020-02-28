use Getopt::Long;
our $VERSION = "1.13.1";

my ($help,$parameter);
GetOptions('-help|h'=>\$help,
			'-p=s'=>\$parameter,
		);

if(!defined($parameter))
{
	help();
}

print "\n  JUMPaq initialized\n\n";
		
my $params=parse_param($parameter);

open(OUTPUT,">$params->{output}") || die "please input the output file\n";

my $db = $params->{database};
if(!defined($db))
{
	print "please input the database file in the parameter file\n";
	exit;
}


#my $db = $ARGV[0];
#my $Report = $ARGV[1];

my $Report = $params->{report_file};

open(REP,$Report) || die "can not open the file: $Report";

my $dbhash = create_dbHash($db);
print "\n";

#$params->{'min_peptide_length'} = 6;
#$params->{'enzyme'} = 'trypsin';
#$params->{'max_peptide_length'} = 30;
#$params->{'max_peptide_hydro'} = 17;
#$params->{'min_peptide_hydro'} = -24;

my %enzymes = (trypsin => "[KR\-](?!P)", chymotrypsin => "[WYF\-](?!P)",
               gluc_nahpo => "[ED\-](?!P)", gluc => "[E\-](?!P)",
               lysc => "[K\-](?!P)", argc => "[R\-](?!P)", aspc => "[D]");


my %phobicity =(A => 16, B => 114.5962, C => 250, D => -249, E => -150, F => 500,
                G => -331, H => -463, I => 441, K => -500, L => 476, M => 332,
                N => -379, O => 114.1472, P => -492, Q => -276, R => -277, S => -285,
                T => -108, V => 302, W => 488, X => 113.1594, Y => 200, Z => 128.6231);

				
my $number_no_theo_peptide = 0;
my $found_in_database = 0;
my $not_found_in_database = 0;

my $summary_line = <REP>;
my $header = <REP>;			
while(<REP>)
{
	my @lines = split(/\t/,$_);
	print "\r  Calculating emPAI for protein: $lines[1]";
	if(!defined($dbhash->{$lines[1]}))
	{
	
		$not_found_in_database++;
		next;
	}
	my $good = calc_emPAI($dbhash->{$lines[1]}->{'sequence'});
	if ($good==0)
	{
		$good = 1;
		next; 
	}
	$found_in_database++;
#	my $PAI = $lines[3]/$good;
	my $PAI = $lines[4]/$good;

#	my $abundance = (10**$PAI) - 1;	
	my $abundance = $PAI * 50000;	

	my @data = split(/\|/,$lines[1]);
	print OUTPUT $data[1],"\t",sprintf("%.0f",$abundance),"\n";
}
print "\n";

print "  $found_in_database proteins were found in the database and were calculated its emPAI value\n";
print "  $not_found_in_database proteins can not be found in the database\n";
 
print "\n  JUMPaq finished\n\n";

				
sub calc_emPAI{
    my ($proteinseq) = @_;

    my $ftpeps = get_peptides_emPAI($proteinseq);
    my $good = 0;
    for my $peptide (@$ftpeps)
	{
        $peptide =~ s/[\-A-Z]+\.([A-Z]+)\.[A-Z\-]+\Z/$1/;
        my @temp = split("", $peptide);
        my $pho = get_Pho(\@temp);
        $good++ if (($pho <= $params->{'max_peptide_hydro'}) && ($pho >= $params->{'min_peptide_hydro'}));
    }
    #$$proteinhash{$protein}{'emPAI_ftpeps'} = $good;
	return $good;
}

sub get_peptides_emPAI{
  my ($sequence) = @_;
  $min_len= $params->{'min_peptide_length'};
  $enzyme = $params->{'enzyme'};
  $max_len = $params->{'max_peptide_length'};
  
  $max_len = 1000 if (!defined($max_len));
  $sequence = "-".$sequence."-";
  $enzyme = lc $enzyme;
  my $pattern = $enzymes{$enzyme} || die "$enzyme is not in Enzyme Hash\n";;
  my @parts = split(/($pattern)/, $sequence);
  my ($i, $arraysize) = (2, scalar(@parts));
  my @tryp_array;
  while ($i < $arraysize){
    if ((length($parts[$i])>=$min_len-1) && (length($parts[$i])<$max_len)){
        my $Nterm = get_Nterm(\@parts, $i);
        my $Cterm = get_Cterm(\@parts, $i);
        push (@tryp_array, "$Nterm$parts[$i]$Cterm");
    }
    $i += 2;
  }
  return \@tryp_array;
}

sub get_Nterm{  #return 2 amino acids on N terminal: both are preceding the internal sequence
    my ($peps, $i) = @_;
    my $Nterm = "";
    my $offset = 1;


    while(length($Nterm) < 3)
	{
        if (defined($$peps[$i-$offset]))
		{
                $Nterm = $$peps[$i-$offset].$Nterm;
        }
        last if $$peps[$i-$offset] =~ /-/;
        $offset++;
    }
	my $length = length($Nterm);
	if ($length>=3)
	{
		$Nterm =~ s/.*(\w\w)\Z/$1\./;
    } elsif ($length==2)
	{
        $Nterm =~ s/([\-\w]{2})/$1\./;
    } elsif ($length==1)
	{
        $Nterm .= "\.";
    }

    return $Nterm;
}
sub get_Cterm{  #returns 2 amino acids on C terminal:  one w/n the internal sequence and one trailing
    my ($peps, $i) = @_;
    my $Cterm = "";
    my $offset = 1;
    while(length($Cterm) < 4)
	{
        if (defined($$peps[$i+$offset]))
		{
            $Cterm .= $$peps[$i+$offset];
        }
        if (!defined($$peps[$i+$offset]))
		{
			for (@$peps)
			{
				print "test $_\n";
            }
            exit;
        }
        last if ($$peps[$i+$offset] =~ /\-/);
        $offset++;
    }
    my $length = length($Cterm);
    if ($length>=4){
            $Cterm =~ s/(\w)([\w\-]{2}).*/$1\.$2/;
    } elsif ($length==3) {
            $Cterm =~ s/(\w)(.*)/$1\.$2/;
    } elsif ($length==2) {
            $Cterm =~ s/(\w)([\w\-])/$1\.$2/;
    } else {
            $Cterm = "\.".$Cterm;
    }
                                                                                              
    return $Cterm;
}

sub get_Pho{  #returns hydrophobicity of given sequence
        shift @_;
  my ($aacids) = @_;
  my $Pho = 0;


  for (@$aacids){
    $Pho += $phobicity{$_} if defined($phobicity{$_});
  }
  return $Pho/100;
}

sub create_dbHash{  #returns a hash that contains name, annotation, and sequence for all proteins in given database

    my ($db) = @_;
	open (DB, "<$db") || die "Cannot open database: $db $!\n";
	my %dbHash;
	my $inProt = 0;
	my ($name, $annotation, $sequence);
	my $number = 0;
	
	my $number_entry = 0;
	while (<DB>){
		chomp($_);
		if (/^>/){
			$number_entry++;
			$inProt = 1;
			$_=~s/\#\#//g;
			$_ =~ s/^>([a-zA-Z0-9\.\_\-\|\:]+)[\s\,\;]*(.*)//;
			$name = $1;
			print "\r  Creating protein database hash: $name";			
			$annotation = $2;
			if (!defined($annotation))
			{
				print $_;
				$annotation = "";
			}
			$dbHash{$name}{'protein'} = $name;
			$dbHash{$name}{'annotation'} = $annotation;
			my $species;
			$annotation =~ s/\s+\[MASS\=\d+\].*\Z//o;
			$annotation =~ s/\s+\Z//;
			if ($annotation =~ /\[[\w\.\-\s]+\]\Z/){
					$annotation =~ s/\[([\w\.\-\s]+)\]\Z//o;
					$species = $1;
			} else {
					$species = "Not specified";
			}
			$dbHash{$name}{'species'} = $species;
			$dbHash{$name}{'number'} = $number;
			$number++;
		} elsif ($inProt){
			$sequence = $_;
			$sequence =~ s/\s//g;
			$sequence =~ s/[\*\#\@\_]//g;
			next if ($sequence eq "");
			$dbHash{$name}{'sequence'} .= $sequence;
		}
	}
	print "\n  $number_entry proteins in the database";
	print "\n";
	return (\%dbHash);
}

sub parse_param {
  my($path) = shift;

	print "  Parsing parameter file\n";
  if(open(P,"< $path")){
    my($line);
    my $phash = {};
    
    while($line = <P>){
      
      my $linehash = {};
      my $comments = "";
      
      if( $line =~ s/\s*([;\#].*)$// ) {$comments = $1;}
      
      $linehash->{Comments} = $comments;
      
      if($line =~ /^(.+?)\s*=\s*(.+)$/){
        my ($key,$data) = ($1,$2);
        $data =~ s/\s+$//o;
#        $linehash->{data} = $data;
        $phash->{$key} = $data;
      }      
   }
    
    close P;
    $self->{PARAMETERS} = $phash;
	return $phash;
  }
  return 0;
  
}

sub help {
		print "\n";
		print "     Usage: JUMPaq -p JUMPaq.params \n";
		print "\n";
		exit;
}

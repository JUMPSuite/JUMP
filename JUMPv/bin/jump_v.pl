#!/bin/env perl

use Getopt::Long;
our $VERSION = "1.13.002";

my ($help,$parameter);
GetOptions('-help|h'=>\$help,
			'-p=s'=>\$parameter,
		);

if(!defined($parameter))
{
	help();
}

print "\nJUMPv initialized\n";
print "Reading parameter file\n";		
my $params=parse_param($parameter);

open(OUTPUT,">$params->{output}") || die "please input the output file\n";

my $file = $params->{quan_file};
if(!defined($file))
{
	print "please input the quan_file in the parameter file\n";
	exit;
}
print "Reading quantification raw file\n";
my ($headline,$peptidehash,$proteinhash) = ReadRawfile($file);
print OUTPUT $headline,"\n";
print "Extracting the peptides/proteins\n";
Extract_data($params,$peptidehash,$proteinhash);

print "JUMPv finished\n\n";

sub ReadRawfile
{
########### Read input file ######################
	my ($file) = shift;
	
	open(IN, '<', $file) || die $!;
	my $database = readline IN; # database line;
	my $headline = readline IN;  # skip first line
	chop $headline;
	my @head_array = split(/\;/,$headline);

	my (%outfilehash, %peptidehash, %proteinhash);

	while(<IN>)
	{	
		chomp;
		my @data = split(/\;/,$_);
		next if ($#data<2);
		my @data1 = split(/\./,$data[0]);
		my $pep = $data1[1];
		print "\rReading peptide: $pep";
		my $protein = $data[1];
		my $out = $data[2];

		$proteinhash{$protein}{$out}=$_;
		$peptidehash{$pep}{$protein}{$out} = $_;		
	}
	close IN;
	print "\n";
	return ($headline,\%peptidehash,\%proteinhash);
}


sub parse_param {
  my($path) = shift;


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

sub Extract_data
{
	my ($param,$peptidehash,$proteinhash) = @_;
	foreach (keys %$param)
	{
		if($_ =~/peptide_/)
		{
			foreach my $peptide (keys %$peptidehash)
			{
				my $peptide_orig = $peptide;
				$peptide =~ s/[^a-zA-Z0-9 ]//g;
				if($peptide eq $param->{$_})
				{
					print "peptide found: $peptide_orig found\n";
					foreach my $protein (keys %{$peptidehash->{$peptide_orig}})
					{					
						print "protein found: $protein found\n";
						
						foreach my $out (keys %{$peptidehash->{$peptide_orig}->{$protein}})
						{
							print OUTPUT $peptidehash->{$peptide_orig}->{$protein}->{$out},"\n";
						}
					}						
				}
			}			
		}
		if($_ =~/protein_/)
		{
			foreach my $protein (keys %$proteinhash)
			{
				if($protein eq $param->{$_})
				{
					print "protein found: $protein\n";				
					foreach my $out (keys %{$proteinhash->{$protein}})
					{				
						print OUTPUT $proteinhash->{$protein}->{$out},"\n";
					}	
				}
			}	
		}
		
	}

}

sub help {
		print "\n";
		print "     Usage: jump -p jump_v.params \n";
		print "\n";
		exit;
}
		

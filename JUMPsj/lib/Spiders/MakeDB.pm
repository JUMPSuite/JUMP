#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::MakeDB
package Spiders::MakeDB; 


######### Job #################################################
#                                                             #
#       **************************************************    #  
#       **** JUMP program           		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2012 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = qw(get_dbfile get_enzyme make_db get_peptides get_Nterm get_Cterm create_dbHash);
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _db_file => $arg{'db_file'},
    };
    bless $self, $class;
    return $self;
}

sub get_dbfile
{
	my $self=shift;
	return $self->{'_db_file'};
}

sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}


sub get_enzyme
{
	my $self=shift;
	
	my %enzymes = (trypsin => "[KR\-](?!P)", 
				   chymotrypsin => "[WYF\-](?!P)",
				   gluc_nahpo => "[ED\-](?!P)",
   				   gluc => "[E\-](?!P)",
                   lysc => "[K\-](?!P)",
 				   argc => "[R\-](?!P)",
				   aspc => "[D]");
				   
	return \%enzymes;			   
}

sub make_db
{
	my $self=shift;
	
	my $db = $self->get_dbfile();
	my $param = $self->get_parameter();
	my $parsedSqs = [];
		
	my ($dbhash, $proteinnum) = $self->create_dbHash($db);
	my $protein_num = scalar keys (%$dbhash);
	print "There are $protein_num in the database\n";


	while (my ($protein, $hash) = each %$dbhash)
	{
		my $head;
		$head->{'id'} = $protein;
		$head->{'desc'} = $$hash{'annotation'};
		$head->{'seq'} = $$hash{'sequence'};

		push @{$parsedSqs}, $head;
		
	}
	if($param->{'inclusion_decoy'}==1)
	{			
		while (my ($protein, $hash) = each %$dbhash)
		{
			my $decoyKey = $protein;
			$decoyKey =~ s/>(\w+)/>##Decoy__$1/;
			my $head;
			$head->{'id'} = $decoyKey;
			$head->{'desc'} = $$hash{'annotation'};
			$head->{'seq'} = reverse($$hash{'sequence'});

			push @{$parsedSqs}, $head;
			
		}		
	}	

	
	return $parsedSqs;
}

sub get_peptides{
    my ($self, $sequence, $min_len, $enzyme) = @_;
	my $enzymes = $self->get_enzyme();
	
    $sequence = "-".$sequence."-";
######### Version1.20 ################
	
 #   $enzyme = lc $enzyme;
 #   my $pattern = $$enzymes{$enzyme} || die "$enzyme is not in Enzyme Hash\n";
	my $pattern = "";
	if($enzyme =~/(\w+)\s+(\w+)\s+(\w+)/)
	{
		$pattern = "[$2\-](?!$3)";
	}
	else
	{
		print "  Please specify a right format of enzyme_info.\n";
	}
	
	my $pattern1 = $pattern;
	$pattern1 =~ s/\-//g;
	$sequence=~s/(\w)($pattern1)/$2$1/g;
    my @parts = split(/($pattern)/, $sequence);
    my ($i, $arraysize) = (2, scalar(@parts));
    my @tryp_array;
    while ($i < $arraysize){
            if (length($parts[$i])>$min_len){
                            my $Nterm = get_Nterm(\@parts, $i);
                            my $Cterm = get_Cterm(\@parts, $i);
                            push (@tryp_array, "$Nterm$parts[$i]$Cterm");
            }
            $i += 2;
    }
    return ($sequence,\@tryp_array);
}

sub get_Nterm{  #return 2 amino acids on N terminal: both are preceding the internal sequence
    my ($peps, $i) = @_;
    my $Nterm = "";
    my $offset = 1;

    while(length($Nterm) < 3){
		if (defined($$peps[$i-$offset])){
			$Nterm = $$peps[$i-$offset].$Nterm;
        }
		last if $$peps[$i-$offset] =~ /-/;
		$offset++;
    }
    my $length = length($Nterm);
    if ($length>=3){
		$Nterm =~ s/.*(\w\w)\Z/$1\./;
    } elsif ($length==2){
		$Nterm =~ s/([\-\w]{2})/$1\./;
	} elsif ($length==1){
		$Nterm .= "\.";
    }
    return $Nterm;
}

sub get_Cterm{  #returns 2 amino acids on C terminal:  one w/n the internal sequence and one trailing
    my ($peps, $i) = @_;
    my $Cterm = "";
    my $offset = 1;

    while(length($Cterm) < 4){
		if (defined($$peps[$i+$offset])){
             $Cterm .= $$peps[$i+$offset];

		}
		if (!defined($$peps[$i+$offset])){
             for (@$peps){
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

sub create_dbHash{  #returns a hash that contains name, annotation, and sequence for all proteins in given database
    my ($self,$db) =  @_;

	open (DB, "<$db") || die "Cannot open database: $db $!\n";

	my %dbHash;
	my $inProt = 0;
	my ($name, $annotation, $sequence);
	my $number = 0;
	while (<DB>)
	{
		chomp($_);
		if (/^>/){
			my $decoy=0;
			$inProt = 1;
			if($_=~/\#\#/)
			{
				$_=~s/\#\#//g;
				$decoy=1;
			}
			$_ =~ s/^>([a-zA-Z0-9\.\_\-\|\:]+)[\s\,\;]*(.*)//;
			if($decoy)
			{
				$name = "##" . $1;
				print $name,"\n";
			}
			else
			{
				$name = $1;
			}
			$annotation = $2;
			if (!defined($annotation)){
				$annotation = "";
			 }
			$dbHash{$name}{'protein'} =  $name;
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
		}
		elsif ($inProt){
			$sequence = $_;
			$sequence =~ s/\s//g;
			$sequence =~ s/[\*\#\@\_]//g;
			next if ($sequence eq "");
			$dbHash{$name}{'sequence'} .= $sequence;
		}
	}
	return (\%dbHash, $number);
}

1;

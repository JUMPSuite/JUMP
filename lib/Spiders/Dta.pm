#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Dta

######### Simulation ##########################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
package Spiders::Dta;

use strict;
use warnings;
use File::Spec;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();


sub new{
    my $class = shift;
    my $dta_file;
    my $instream;
    my $charge;
    my $precMZ;
    my $mz_int_hash = {};
    my $mz_array;
    my $int_array;   
    my $self;
    if(scalar(@_) > 1) {
	$dta_file = shift;    
	$charge = shift;	     
	$precMZ = shift;	     
	$mz_array = shift;    
	$int_array = shift;   

	unless( scalar(@$mz_array) == scalar(@$int_array) ) {
	    die( "mz_array and int_array must have the same length" );
	}

	for( my $i = 0; $i < scalar(@$mz_array); ++$i ) {
	    $mz_int_hash->{$$mz_array[$i]} = $$int_array[$i];
	}

	$self = {
	    '_dta_file'  => $dta_file,
	    '_charge' => $charge,	     
	    '_precMZ' => $precMZ,	     
	    '_mz_int_hash' => $mz_int_hash, 
	    '_mz_array' => $mz_array,    
	    '_int_array' => $int_array
	};
    }
    elsif(scalar(@_) == 1){
	$dta_file = shift;
	$self = process_dtafile($dta_file);
	my ($vol,$dir,$file) = File::Spec->splitpath($dta_file);
	$self->{'_dta_file'} = $file;
    }

    bless $self, $class;

    $self->parse_dtafile_name();
    return $self;
}

# sub set_dta_file
# {
# 	my ($self,$dtafile) = @_;
# 	return $self->{'_dta_file'} = $dtafile;	
# }

sub get_dta_file
{
    my $self = shift;
    return $self->{'_dta_file'};
}

sub parse_dtafile_name
{
    my $self=shift;
    my $dtafile = $self->get_dta_file();
    $dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
    my ($scan, $specscan, $specscan2, $charge) = ($3, $3, $4, $5); 
    $self->{'_scan'} = $scan;
    $self->{'_charge'} = $charge;
}

sub get_scan_num
{
    my $self = shift;
    return $self->{'_scan'};
}

sub get_charge
{
    my $self = shift;
    return $self->{'_charge'};
}

sub process_dtafile
{
    my $dtafile = shift;
    my $self = {};

# define the mz hash for storing the mz and intensity
    my %mz_hash;
# open each dta file to load the data 
    open(DTAFILE,$dtafile) || die "can not open the dta file: $dtafile";
# get the precusor mz and charge 
    my $prec_mz_charge = <DTAFILE>;
    chomp $prec_mz_charge;
    my ($prec_mz,$prec_charge) = split(/\s+/,$prec_mz_charge);
    my @mz_array;
    my @int_array;
    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
	my @data =split(/\s+/,$_);
	$mz_hash{$data[0]} = $data[1];
	push (@mz_array,$data[0]);
	push (@int_array,$data[1]);
    }
    close(DTAFILE);
    
    $self->{'_charge'} = $prec_charge;
    $self->{'_precMZ'} = $prec_mz;
    $self->{'_mz_int_hash'} = \%mz_hash;	$
	self->{'_mz_array'}=\@mz_array;
    $self->{'_int_array'} = \@int_array;
    return $self
}

sub to_string_hashorder
{
    my ($self,$precision) = @_;
    my $mz_int_hash = $self->get_mz_int_hash();
    my $repr = $self->get_prec_mz() . " " . $self->get_charge() . "\n";

    if( defined($precision) ) {
	my $fmt = "%" . $precision . "f";
	foreach my $mz (sort {$a<=>$b} keys %{$mz_int_hash}) {
	    my $intensity = sprintf($fmt,$mz_int_hash->{$mz});
	    my $fmz = sprintf($fmt,$mz);
	    $repr .= $fmz . " " . $intensity . "\n";
	}
    }
    else {
	foreach my $mz (sort {$a<=>$b} keys %{$mz_int_hash}) {
	    $repr .= $mz . " " . $mz_int_hash->{$mz} . "\n";
	}
    }

    return $repr;
}

sub to_string_arrorder
{
    my ($self,$precision) = @_;
    my $int_array = $self->get_int_array();
    my $mz_array = $self->get_mz_array();
    my $repr = $self->get_prec_mz() . " " . $self->get_charge() . "\n";

    if( defined($precision) ) {
	my $fmt = "%" . $precision . "f";
	for( my $i = 0; $i <= $#$mz_array; ++$i ) { 
	    my $intensity = sprintf($fmt,$$int_array[$i]);
	    my $fmz = sprintf($fmt,$$mz_array[$i]);
	    $repr .= $fmz . " " . $intensity . "\n";
	}
    }
    else {
	for( my $i = 0; $i <= $#$mz_array; ++$i ) { 
	    $repr .= $$mz_array[$i] . " " . $$int_array[$i] . "\n";
	}
    }

    return $repr;
}

sub write_dta_file_hashorder
{
    my ($self,$precision) = @_;
    my $dtafile = $self->get_dta_file();
    
    # if(-e $dtafile)
    # {
    # 	system(qq(rm $dtafile));
    # }
    open (DTA, ">$dtafile") || die "can not open the dta file: $dtafile\n";
    print DTA $self->to_string_hashorder();
    close(DTA);
}

sub write_dta_file_arrorder
{
    my ($self,$precision) = @_;
    my $dtafile = $self->get_dta_file();
    
    open (DTA, ">$dtafile") || die "can not open the dta file: $dtafile\n";
    print DTA $self->to_string_arrorder();
    close(DTA);
}

# sub set_prec_mz
# {
# 	my ($self,$prec_mz) = @_;
# 	$self->{'_precMZ'} = $prec_mz;
# }

# sub set_charge
# {
# 	my ($self,$charge) = @_;
# 	$self->{'_charge'} = $charge;
# }

# sub set_mz_int_hash
# {
# 	my ($self,$mz_int_hash) = @_;
# 	$self->{'_mz_int_hash'} = $mz_int_hash;
# }

sub get_prec_mz
{
    my $self = shift;
    return $self->{'_precMZ'};
}

sub get_mz_int_hash
{
    my $self = shift;
    return $self->{'_mz_int_hash'};
}

sub get_mz_array
{
    my $self = shift;
    return $self->{'_mz_array'};
}

sub get_int_array
{
    my $self = shift;
    return $self->{'_int_array'};
}

sub update 
{
    my ($self,$new_mz,$new_int) = @_;
    $self->{'_int_array'} = $new_int;
    $self->{'_mz_array'} = $new_mz;
    $self->{'_mz_int_hash'} = {};
    for( my $i = 0; $i < scalar(@$new_mz); ++$i ) {
	$self->{'_mz_int_hash'}->{$$new_mz[$i]} = $$new_int[$i];
    }
}

1;

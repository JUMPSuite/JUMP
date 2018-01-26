#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Dtas

package Spiders::Dtas;

use strict;
use warnings;

sub new
{
    my ($class,$backend) = @_;
    my $self = { 'backend' => $backend };
    bless ($self,$class);
    return $self;
}

#-------------------------------------------------------------------------------------------
1;

sub add_dta 
{
    my ($self,$dta)=@_;
    $self->{'backend'}->store_dta( $dta->get_dta_file(), $dta );
}

sub get_dta
{
    my ($self,$dta_file) = @_;
    return $self->{'backend'}->load_dta( $dta_file );
}

sub list_dta 
{
    my ($self) = @_;
    return $self->{'backend'}->keys();
}

sub print_dtas 
{
    my ($self,$outfilename,$dtasfiles)=@_;
    open( my $outf, ">>$outfilename" );
    foreach my $dta_file (@$dtasfiles) 
    {
	my $dta = $self->get_dta($dta_file);
	print $outf "$dta_file ",$dta->get_prec_mz()," ",$dta->get_charge(),"\n";
	print $outf join( " ", @{$dta->get_mz_array()} ),"\n";
	print $outf join( " ", @{$dta->get_int_array()} ),"\n";
    }
    close( $outf );
}


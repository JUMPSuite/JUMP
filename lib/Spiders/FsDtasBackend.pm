## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::FsDtasBackend

# A class to be a backend for storing Dtas data on the filesystem.
# Files are stored in the directory pointed to by $path

package Spiders::FsDtasBackend;

use strict;
use warnings;
use File::Spec;
use Spiders::Dta;
use Carp;

sub new {
    my ($class,$path) = @_;
    my $fh = {};
    foreach my $f (glob(File::Spec->join($path,'*.dta'))) {
	my ($vol,$path,$base) = File::Spec->splitpath($f);
	$fh->{$base} = 1;
    }

    my $self = {'path' => $path,
		'dta_files' => $fh};
    bless $self, $class;
    return $self;
}

sub load_dta {
    my ($self,$dta_file) = @_;
    my $abspath = File::Spec->join($self->{'path'},$dta_file);
    unless(defined($self->{'dta_files'}->{$dta_file})) {
	croak( "DTA file $dta_file not found in $self->{'path'}\n" );
    }

    return Spiders::Dta->new($abspath);
}

sub store_dta {
    my ($self,$dta_file,$dta_obj) = @_;
    my $abspath = File::Spec->join($self->{'path'},$dta_file);
    $self->{'dta_files'}->{$dta_file} = 1;
    open(my $outf, ">$abspath");
    print $outf $dta_obj->to_string_arrorder();
    close( $outf );
}

sub keys {
    my ($self) = @_;
    my $href = $self->{'dta_files'};
    my @retval = keys(%$href);
    return \@retval;
}

1;

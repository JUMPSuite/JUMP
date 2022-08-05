## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::IdxDtasBackend

# A class to be a backend for storing Dtas data on the filesystem.
# Files are frozen and stored in one big flat file

package Spiders::IdxDtasBackend;

use strict;
use warnings;
use Exporter;
use File::Spec;
use Spiders::Dta;
use Storable qw(store retrieve freeze thaw);
use Carp;

sub new {
    my ($class,$path,$mode) = @_;
    my $fh = {};
    my $idxfilename = File::Spec->join($path,"dta.idx");
    my $mapfilename = File::Spec->join($path,"dta.map");
    my $mapfile;
    my $idxfile;
    if( -e $idxfilename && -e $mapfilename ) {
	$mapfile = retrieve($mapfilename);
    }
    else {
	$mapfile = {};
	open( $idxfile, "+>", $idxfilename );
	binmode( $idxfile );
    }

    open( $idxfile, "+<", $idxfilename );
    binmode( $idxfile );
    my $self = {'path' => $path,
		'idxfile' => $idxfile,
		'mapfile' => $mapfile,
		'mode' => $mode};
    bless $self, $class;
    return $self;
}

sub DESTROY {
    my $self = shift;
    if( $self->{'mode'} eq 'create' ) {
	store( $self->{'mapfile'}, File::Spec->join($self->{'path'},'dta.map') );
    }
}

sub load_dta {
    my ($self,$dta_file) = @_;
    if( !defined($self->{'mapfile'}->{$dta_file}) ) {
	croak( "DTA file $dta_file not found in IDX file" );
    }
    seek( $self->{'idxfile'}, $self->{'mapfile'}->{$dta_file}, 0 );
    my $bytes;
    read( $self->{'idxfile'}, $bytes, length(pack('L',0)) );
    my @nbytes = unpack( 'L', $bytes );
    read( $self->{'idxfile'}, $bytes, $nbytes[0] );
    return thaw($bytes);
}

sub store_dta {
    my ($self,$dta_file,$dta_obj) = @_;
    my $bytes = freeze($dta_obj);
    my $size = length($bytes);
    if( defined($self->{'mapfile'}->{$dta_file}) ) {
	my $result;
	seek( $self->{'idxfile'}, $self->{'mapfile'}->{$dta_file}, 0 );
	read( $self->{'idxfile'}, $result, length(pack('L',0)) );
	$result = unpack( 'L', $result );
	unless( $result >= $size ) {
	    croak( "cannot replace $dta_file with a larger entry: $size > $result" );
	}
	seek( $self->{'idxfile'}, $self->{'mapfile'}->{$dta_file}, 0 );
	print {$self->{'idxfile'}} pack('L',$size);
	print {$self->{'idxfile'}} $bytes;
    }
    else {
	$self->{'mapfile'}->{$dta_file} = tell($self->{'idxfile'});
	print {$self->{'idxfile'}} pack('L',$size);
	print {$self->{'idxfile'}} $bytes;	
    }
}

sub keys {
    my ($self) = @_;
    my $href = $self->{'mapfile'};
    my @retval = keys(%$href);
    return \@retval;
}

1;

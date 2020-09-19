## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::FsDtasBackend

# A class to be a backend for storing Dtas data on the filesystem.
# Files are stored in the directory pointed to by $path

package Spiders::LmdbDtasBackend;

use strict;
use warnings;
use File::Spec;
use Spiders::Dta;
use LMDB_File qw(:flags :cursor_op);

sub new {
    my ($class,$path,$mode) = @_;
    my $env = LMDB::Env->new($path,
			     {mapsize => 128*1024*1024*1024,
			      maxdbs => 20,
			      mode => 600});
    my $transaction = $env->BeginTxn();
    my $db = $transaction->OpenDB({dbname => "jump_dtas_db",
				   
    my $self = {'path' => $path,
		
    bless $self, $class;
    return $self;
}

sub load {
    my ($self,$dta_file) = @_;
    my $abspath = File::Spec->join($self->{'path'},$dta_file);
    unless(defined($self->{'dta_files'}->{$abspath})) {
	croak( "DTA file $dta_file not found in $self->{'path'}\n" );
    }

    return Spiders::Dta->new($abspath);
}

sub store {
    my ($self,$dta_file,$dta_obj) = @_;
    my $abspath = File::Spec->join($self->{'path'},$dta_file);
    $self->{'dta_files'}->{$abspath} = 1;
    open(my $outf, ">$abspath");
    print $outf $dta_obj->to_string();
    close( $outf );
}

sub keys {
    my ($self) = @_;
    my $href = $self->{'dta_files'};
    return keys(%$href);
}

1;

######### MakeJobManager ######################################
#                                                             #
#       **************************************************    #  
#       **** Job Manager implementation using GNUmake ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 2019 - Alex Breuer	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****alex.breuer@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::MakeJobManager;

use strict;
use File::Temp;
use Spiders::BatchSystem;

sub new {
    my $class = shift;
    my $options = shift || {};
    my $self = {};
    bless $self,$class;

    $self->{'myDir'} = File::Temp->newdir( TEMPLATE => 'XXXXXXXXXXX',
					   UNLINK => (defined($options->{'DEBUG'})) );
    return $self;
}

sub _generateMakefile {
    my $self = shift;
    my @jobs = @_;
    my $batchSystem = new Spiders::BatchSystem();

    my %rules;
    my $i = 0;
    foreach my $j (@jobs) {
	my $cmd = "@" . $batchSystem->getBatchCmd( $j->{'toolType'} ) . ' "' . $j->{'cmd'} . '"';
	my $sfile = File::Temp->new( DIR => $self->{'myDir'} );
	open( my $sout, ">".$sfile->filename );
	print $sout "#!/bin/bash\n$cmd &> /dev/null\n";
	close( $sout );
	chmod( 0700, $sfile->filename );
	$rules{"$i-cmd-run"} = $sfile;
	$i += 1;
    }

    my $tfile = File::Temp->new( TEMPLATE => 'XXXXXX',
				 DIR => $self->{'myDir'} );
    open( my $outf, ">".$tfile->filename );
    print $outf 'all: ' . join( ' ', keys(%rules) ) . "\n";
    foreach my $k (keys(%rules)) {
	print $outf "$k:\n\t" . $rules{$k} ." ; echo finished job" . ($k =~ /\d+/) . "\n";
    }
    close( $outf );
    return $tfile;
}

sub runJobs {
    my $self = shift;
    my @jobs = @_;

    my $mfile = $self->_generateMakefile( @jobs );
    system( "make -j -f " . $mfile->filename );
}

1;

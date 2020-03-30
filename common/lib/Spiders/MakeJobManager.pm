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
use Carp;
use List::Util qw[min];
use File::Temp;
use Spiders::BatchSystem;
use Spiders::Config;
use Spiders::ClusterConfig;

sub new {
    my $class = shift;
    my $params = shift;
    my $options = shift || {};
    my $self = {};
    bless $self,$class;

    $self->{'config'} = new Spiders::Config();
    $self->{'params'} = $params;
    $self->{'unroll'} = $self->{'config'}->get( 'batch_job_unroll' );
    $self->{'debug'} = defined($options->{'DEBUG'});
    $self->{'myDir'} = File::Temp->newdir( TEMPLATE => ($self->{'debug'} ?
								'MAKE_JOB_MANAGER_TEMPDIRXXXX' : 'MakeJobManager.XXXXXXXXXXX'),
					   UNLINK => (!$self->{'debug'}),
					   DIR => '.' );
    return $self;
}

sub _generateSMPMakefile {
    my $self = shift;
    my @jobs = @_;
    
    my %rules;
    my $njobs = scalar(@jobs);
    for( my $i = 0; $i < $njobs; $i += 1 ) {
	my $cmd = $jobs[$i]->{'cmd'};
	chomp($cmd);
	$cmd =~ s/\$/\$\$/g;
	$rules{"$i-cmd-run"} = "( $cmd ) > /dev/null 2>&1  && echo finished job " . $i;
    }

    my $tfile = File::Temp->new( TEMPLATE => 'SMPXXXXXX',
				 DIR => $self->{'myDir'},
				 UNLINK => 0 );
    open( my $outf, ">".$tfile->filename );
    print $outf 'all: ' . join( ' ', keys(%rules) ) . "\n";
    foreach my $k (keys(%rules)) {
	print $outf "$k:\n\t@" .  $rules{$k} . "\n";
    }
    close( $outf );
    return ($tfile,\%rules);

}

sub _generateBatchMakefile {
    my $self = shift;
    my @jobs = @_;
    my $batchSystem = new Spiders::BatchSystem();

    my %cmds;
    my $njobs = scalar(@jobs);


    for( my $i = 0; $i < $njobs; $i += $self->{'unroll'} ) {
	my $cmd = $jobs[$i]->{'cmd'};
	chomp($cmd);
	$cmd =~ s/\$/\$\$/g;
	(my $sfile, my $srules) = $self->_generateSMPMakefile( @jobs[$i..min(scalar(@jobs)-1,
									     ($i + $self->{'unroll'}))] );
	$cmds{$i} = '('.$batchSystem->getBatchCmd( $jobs[$i]->{'toolType'} ) . ' "make -j 1 -f ' . $sfile .'" &> /dev/null  && echo finished job ' . $i . ')';
    }

    my $config = new Spiders::Config();
    my $lag = $config->get('batch_dispatch_lag');
    my $tfile = File::Temp->new( TEMPLATE => 'XXXXXX',
				 DIR => $self->{'myDir'} );
    open( my $outf, ">".$tfile->filename );
    print $outf "all:\n";
    foreach my $k (keys(%cmds)) {
	print $outf "\t\@sleep $lag ; $cmds{$k} &\n";
    }
    print $outf "\t\@wait\n";
    close( $outf );
    return ($tfile,\%cmds);
}

sub runJobs {
    my $self = shift;
    my @jobs = @_;
    my $clusterConfig = Spiders::ClusterConfig::getClusterConfig($self->{'config'},$self->{'params'});

    my $code = 0;
    my $mfile;
    my $handle;
    print "  Submitting " . scalar(@jobs) . " jobs for search\n";
    if( $clusterConfig eq Spiders::ClusterConfig::SMP ) {
	($mfile,my $rules) = $self->_generateSMPMakefile( @jobs );
	open( $handle, "make -f " . $mfile->filename . ' -j ' . $self->{'config'}->get( 'max_search_worker_procs' ) . '|' );
    }
    elsif( $clusterConfig eq Spiders::ClusterConfig::CLUSTER ) {
	($mfile,my $rules) = $self->_generateBatchMakefile( @jobs );	
	open( $handle, "make -f " . $mfile->filename . ' -j ' . $self->{'config'}->get( 'max_dispatch_worker_procs' )  . '|');
    }

    $| = 1;
    print "\r  0 jobs are done";
    for( my $jobsDone = 1; <$handle>; $jobsDone += 1 ) {
	print "\r  $jobsDone jobs are done";
    }
    print "\n";
    $| = 0;

    close $handle;
    $code = $?;

    if( $code ne 0 ) {
	$self->{"myDir"}->unlink_on_destroy( 0 );
	$mfile->unlink_on_destroy( 0 );
	croak( "make retured with nonzero exit code; see makefiles in " . $mfile . "\n" );
    }
}

1;

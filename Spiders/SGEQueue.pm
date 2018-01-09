package Spiders::SGEQueue;

use strict;
use File::Spec;

sub new {
    my ($class,%arg) = @_;
    my $self = {};
    bless $self, $class;
    return $self;
}

sub submit_job {
    my ($path,$name,$cmd) = @_;

    my $file = File::Spec->catfile($path,"$name.sh");
    open (JOB, ">$file") or die "Cannot creat a job file\n";
    print JOB "#!/bin/bash\n";
    print JOB "#\$ -N $name\n";
    print JOB "#\$ -e $path/$name.e\n";
    print JOB "#\$ -o $path/$name.o\n";
    print JOB $cmd;
    close (JOB);
    my $command = qq(qsub -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G $file);
    my $job = lc(qx[$command]);
    chomp($job);
    if ($job =~ /job (\d+)/) {
	$job = $1;
    }
    else {
	warn "could not parse qsub output";
    }
    return $job
}

sub get_running_jobs {
    my ($self,$jobshash) = @_;
    my $command =  "qstat";
    my $jobStatus = qx[$command];
    my @lines = split( /\n/, $jobStatus );
    shift @lines;
    my @retval;
    foreach my $l (@lines) {
	my @toks = split( /\s+/, $l );
	if( defined($jobshash->{$toks[0]}) ) {
	    push( @retval, $toks[0] );
	}
    }
    return \@retval;
}

1;

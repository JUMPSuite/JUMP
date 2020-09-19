package Spiders::LSFQueue;

use strict;
use File::Spec;

sub new {
    my ($class,%arg) = @_;
    my $self = {};
    bless $self, $class;
    return $self;
}

sub submit_job {
    (my $self,my $path,my $name,my $cmd) = @_;

    my $file = File::Spec->catfile($path,"$name.sh");
    open (JOB, ">$file") or die "Cannot create a job file\n";
    print JOB "#!/bin/bash\n";
    print JOB "#BSUB -P prot\n";
    print JOB "#BSUB -q standard\n";
    print JOB "#BSUB -R \"rusage[mem=20000]\"\n";
    print JOB "#BSUB -eo $path/$name.e\n";
    print JOB "#BSUB -oo $path/$name.o\n";
    print JOB $cmd;
    close (JOB);
    my $command = qq(cd $path && bsub < $name.sh);
    my $job = qx[$command];
    chomp($job);
    if ($job =~ /Job \<(\d+)\>/) {
	$job = $1;
    }
    else {
	warn "could not parse bsub output: $job";
    }
    return $job
}

sub get_running_jobs {
    my ($self,$jobshash) = @_;
    my $command =  "bjobs -noheader";
    my $jobStatus = qx[$command];
    my @lines = split( /\n/, $jobStatus );
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

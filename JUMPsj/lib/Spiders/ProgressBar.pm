## Module name: Spiders::ProgressBar

####################### ProgressBar ###########################
#                                                             #
#       **************************************************    #  
#       **** JUMP program                             ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****xusheng.wang@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::ProgressBar;

use strict;

sub new {
    my $class = shift;
    my $ntot = shift;
    my $minorTicks = shift || 5;
    my $majorTicks = shift || 25;
    my $cycleInterval = shift || 50;

    my $self = {"ntot" => $ntot,
		"majorTicks" => $majorTicks,
		"minorTicks" => $minorTicks,
		"curPercent" => 0,
		"counter" => 0,
		"cycleCharArr" => ['-','\\','|','/'],
		"cycleInterval" => $cycleInterval,
		"cycleI" => 0,
		"backupCursor" => 0};
    
    bless($self,$class);
    return $self;
}

sub incr {
    my $self = shift;
    my $amt = shift || 1;
    $self->{"counter"} += $amt;
    if( $self->{"counter"} % $self->{"cycleInterval"} == 0 ) {
	$| = 0;
	if( $self->{"backupCursor"} == 1 ) {
	    print "\b";
	}
	print $self->{"cycleCharArr"}->[$self->{"cycleI"}];
	$self->{"cycleI"} = ($self->{"cycleI"} + 1) % scalar(@{$self->{"cycleCharArr"}});
	$self->{"backupCursor"} = 1;
	$| = 1;
    }
    if( ($self->{"counter"}/$self->{"ntot"})*100 >= $self->{"curPercent"} ) {
	$| = 1;
	if( $self->{"curPercent"} % $self->{"majorTicks"} == 0 ) {
	    if( $self->{"backupCursor"} == 1 ) {
		print "\b";
	    }
	    print $self->{"curPercent"},"%";
	    $self->{"backupCursor"} = 0;
	    if( $self->{"curPercent"} == 100 ) { print "\n"; }
	}
	else {
	    if( $self->{"backupCursor"} == 1 ) {
		print "\b";
	    }
	    print ". ";
	}
	$self->{"curPercent"} += $self->{"minorTicks"};
	$| = 0;
    }
}

1;

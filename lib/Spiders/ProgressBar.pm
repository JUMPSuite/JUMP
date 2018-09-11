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

    my $self = {"ntot" => $ntot,
		"majorTicks" => $majorTicks,
		"minorTicks" => $minorTicks,
		"curPercent" => 0,
		"counter" => 0};
    
    bless($self,$class);
    return $self;
}

sub incr {
    my $self = shift;
    my $amt = shift || 1;
    $self->{"counter"} += $amt;
    if( ($self->{"counter"}/$self->{"ntot"})*100 >= $self->{"curPercent"} ) {
	$| = 1;
	if( $self->{"curPercent"} % $self->{"majorTicks"} == 0 ) {
	    print $self->{"curPercent"},"%";
	    if( $self->{"curPercent"} == 100 ) { print "\n"; }
	}
	else {
	    print ".";
	}
	$self->{"curPercent"} += $self->{"minorTicks"};
	$| = 0;
    }
}

1;

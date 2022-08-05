#!/bin/env perl

######### Job #################################################
#                                                             #
#       **************************************************    #  
#       **** Dta directory management   		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2018 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::DtaDir;

use strict;
use Carp;
use File::Spec;
use File::Temp;

sub new {
    my ($class,$linkLocation,$behindTheScenesLocation,$keep) = @_;
    my $self = {'linkLocation'=>$linkLocation,
		'behindTheScenesLocation'=>$behindTheScenesLocation,
		'keep'=>$keep};
    if(defined($behindTheScenesLocation)) {
	$self->{'tempDirObj'} = File::Temp->newdir( ".jumpXXXXXXXXXXXXXX",
						    DIR=>$behindTheScenesLocation,
						    UNLINK=>!defined($keep) );
	symlink( $self->{'tempDirObj'}->dirname, $linkLocation );
    }
    else {
	mkdir($linkLocation);
    }
    bless $self,$class;
    return $self;
}

sub DESTROY {
    my $self = shift;
    if(!defined($self->{'keep'}) && defined($self->{'tempDirObj'})) {
	unlink( $self->{'linkLocation'} );
    }
}

1;

######### Which################################################
#                                                             #
#       **************************************************    #  
#       **** Configuration utilities                  ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 2020 - Alex Breuer	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****alex.breuer@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Which;

use strict;
use Carp;

sub which {
    my $cmd = shift;
    my $path = qx[which $cmd];
    if( $? eq 0 ) {
	chomp($path);
	return $path;
    }
    else {
	croak("could not frin $cmd on the path");
    }
}
1;

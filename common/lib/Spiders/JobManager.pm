######### JobManager ##########################################
#                                                             #
#       **************************************************    #  
#       **** Job Manager utilities                    ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 2019 - Alex Breuer	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****alex.breuer@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::JobManager;

use strict;
use Spiders::MakeJobManager;

sub newJobManager {
    my $self = Spiders::MakeJobManager->new({'DEBUG'=>''});
    return $self;
}

1;

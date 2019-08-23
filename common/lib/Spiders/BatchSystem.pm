######### BatchSystem ##########################################
#                                                             #
#       **************************************************    #  
#       **** Batch system utilities                   ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 2019 - Alex Breuer	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****alex.breuer@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::BatchSystem;

use strict;
use Spiders::Config;

use constant {
    JUMP_SEARCH => "jump_search",
    JUMP_DATABASE => "jump_database",
    JUMP_QUANTIFICATION => "jump_quantification",
    JUMP_LOCALIZATION => "jump_localization",
    RUNSEARCH_SHELL => "runsearch_shell"
};

sub new {
    my $class = shift;
    my $self = {};
    bless $self,$class;

    return $self;
}

sub getBatchCmd {
    (my $self, my $tool) = @_;
    my $config = new Spiders::Config();

    my $cmd = $config->get($tool . '_batch_cmd');
    return (defined($cmd) ? $cmd : $config->get('default_batch_cmd'));
}

1;

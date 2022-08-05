######### ClusterConfig #################################################
#                                                             #
#       **************************************************    #  
#       **** Cluster Configuration utilities                  ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 2019 - Alex Breuer	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****alex.breuer@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::ClusterConfig;

use constant {
    CLUSTER => "CLUSTER",
    SMP => "SMP"
};

sub getClusterConfig {
    (my $config, my $params) = @_;
    my $clusterMode;
    if(defined($params->{'cluster'})) {
	warn("param file cluster setting is deprecated and will be removed");
	if(($params->{'cluster'} || $config->get('cluster')) && !($params->{'cluster'} && $config->get('cluster'))) {
	    warn("param and config do not agree on cluster mode; respecting param definition");
	}
	$clusterMode = ($params->{'cluster'} ? CLUSTER : SMP);
    }
    else {
	$clusterMode = ($config->get("cluster") ? CLUSTER : SMP);
    }
    return $clusterMode;
}

1;

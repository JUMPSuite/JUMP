######### Config ##############################################
#                                                             #
#       **************************************************    #  
#       **** Configuration utilities                  ****    #     
#       ****				              ****    #  
#       ****Copyright (C) 2018 - Alex Breuer	      ****    #     
#       ****all rights reserved.	              ****    #  
#       ****alex.breuer@stjude.org	              ****    #  
#       ****				              ****    #  
#       ****				              ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Config;

use strict;
use Storable;
use File::Spec;

use constant {
    SITE_CFG => $ENV{'JUMP_CONFIG_PATH'},
    USER_CFG => File::Spec->join($ENV{'HOME'},'.jump','config')
};

sub new {
    my $class = shift;
    my $self = {};
    bless $self,$class;

    $self->{'userConfig'} = $self->getUserCfg();
    $self->{'globalConfig'} = $self->getGlobalCfg();

    return $self
}

sub getGlobalCfg {
    return retrieve(SITE_CFG);
}

sub getUserCfg {
    my $cfg = {};
    open(my $infile, "<".USER_CFG) || return $cfg;
    while(<$infile>) {
	chomp($_);
	if( length($_) > 0 ) {
	    my @toks = split( /\s+/, $_ );
	    my $k = shift(@toks);
	    $cfg->{$k} = join( ' ', @toks );
	}
    }
    return $cfg;
}

sub keys {
    my $self = shift;
    my %kh;
    foreach my $k (keys %{$self->{'userConfig'}}) {
	$kh{$k} = 1;
    }
    foreach my $k (keys %{$self->{'globalConfig'}}) {
	$kh{$k} = 1;
    }
    return keys %kh;
}

sub provenance {
    my ($self,$key) = @_;
    if(defined($self->{'userConfig'}->{$key})) {
	return USER_CFG;
    }
    elsif(defined($self->{'globalConfig'}->{$key})) {
	return SITE_CFG;
    }
    else {
	return undef;
    }    
}

sub get {
    my ($self,$key) = @_;
    if(defined($self->{'userConfig'}->{$key})) {
	return $self->{'userConfig'}->{$key};
    }
    elsif(defined($self->{'globalConfig'}->{$key})) {
	return $self->{'globalConfig'}->{$key};
    }
    else {
	return undef;
    }
}

sub put {
    my ($self,$key,$value) = @_;
    $self->{'globalConfig'}->{$key} = $value;
    store($self->{'globalConfig'},$ENV{'JUMP_CONFIG_PATH'});
}

1;

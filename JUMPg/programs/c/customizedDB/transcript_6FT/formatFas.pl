#!/bin/env perl 
use lib "$ENV{'JUMP_G_ROOT'}/programs/c";
use PrimarySeq;
$ps=PrimarySeq->new;
$seq=$ps->parseFas($ARGV[0]);
foreach $id (keys %$seq) {
	print ">$id\n$seq->{$id}\n";
}


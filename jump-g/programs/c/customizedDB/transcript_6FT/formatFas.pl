#!/usr/bin/perl -I $JUMP_G_ROOT/c
use PrimarySeq;
$ps=PrimarySeq->new;
$seq=$ps->parseFas($ARGV[0]);
foreach $id (keys %$seq) {
	print ">$id\n$seq->{$id}\n";
}


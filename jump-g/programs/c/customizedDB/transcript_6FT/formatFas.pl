#!/usr/bin/perl -I /home/yli4/development/JUMPg/HPC_version/090118_fullversion/programs/c
use PrimarySeq;
$ps=PrimarySeq->new;
$seq=$ps->parseFas($ARGV[0]);
foreach $id (keys %$seq) {
	print ">$id\n$seq->{$id}\n";
}


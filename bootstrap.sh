#!/bin/bash
show_help() {
    cat <<EOF 
This is the JUMP proteomics pipeline bootstrapping script.  Execute
with no arguments for a standard installation.  An installation of
conda or minicoda is a prerequisite.  

Bootstrapping will create a conda environment in this directory for
use with JUMP.
 
EOF
}

show_success_message() {
    cat <<EOF


JUMP sucessfully bootstrapped.  

For convenience, you may add JUMP to your path:

PATH=\$PATH:$PWD/JUMP/bin
EOF
}

if [[ "$@" =~ .*--help.* ]] ; then
    show_help
    exit 1 
fi

echo "creating conda environment $PWD/conda"
if [ ! -e $(which conda) ] ; then
    echo "conda not found on $PATH ; please ensure conda (or miniconda) is installed and in your PATH variable"
    exit 255;
fi
if [ ! -e $PWD/conda ] ; then
    mkdir $PWD/conda
fi

conda create -y -p $PWD/conda -c conda-forge -c bioconda -c r perl-mime-base64 perl-data-dumper perl-class-std perl-file-spec perl-pod-usage perl-class-std perl-file-spec perl-pod-usage perl-statistics-distributions perl-file-temp perl-excel-writer-xlsx perl-statistics-r perl-file-find perl-http-date perl-math-bigint perl-list-util perl-mime-base64 perl-getopt-long perl-data-dumper perl-parallel-forkmanager perl-clone perl-app-cpanminus r-limma r-fnn r-mass 

if [ $? -ne 0 ] ; then
    echo "Error in conda installation; aborting."
    exit 254
fi

echo "installing cpan modules"
$PWD/conda/bin/cpanm HTTP::Message~"<= 6.20" File::Copy File::Basename Scalar::Util LWP::UserAgent Set::Partition Sys::Hostname

if [ $? -ne 0 ] ; then
    echo "Error in CPAN module installation; aborting."
    exit 253
fi

echo "configuring JUMP"
$PWD/conda/bin/perl Makefile.PL $@ "PERL_BIN=$PWD/conda/bin"
show_success_message

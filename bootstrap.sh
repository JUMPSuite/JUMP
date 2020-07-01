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
add

PATH=$PWD/JUMP/bin:\$PATH

to your .bashrc file or

setenv PATH $PWD/JUMP/bin:\$PATH

to your .cshrc file
EOF
}

while getopts 'dh' OPTION; do
  case "$OPTION" in
      h) show_help
	  exit 1
	  ;;
      d) debug=1
	  ;;
      ?) show_help
	  exit 1
	  ;;
  esac
done
shift "$(($OPTIND -1))"
if [ -n "$debug" ] ; then echo args: $@ ; fi

echo "creating conda environment $PWD/conda"
if [ ! -e $(which conda) ] ; then
    echo "conda not found on $PATH ; please ensure conda (or miniconda) is installed and in your PATH variable"
    exit 255
fi
if [ ! -e $PWD/conda ] ; then
    mkdir $PWD/conda
fi

case $(uname) in 
    Darwin)
	CUDA="" ;;
    Linux)
	CUDA=$(nvcc --version | grep -o 'release ...' | tr -d '[a-z ]') ;;
    *)
	echo "error: not a MacOS or Linux installation"
esac

conda create -p $PWD/conda -y \
  -c bioconda \
  -c conda-forge \
  -c defaults \
  -c pytorch \
  bioconductor-limma=3.40.0 \
  gsl=2.5 \
  harfbuzz=2.4.0 \
  icu=64.2 \
  isl=0.19 \
  jpeg=9c \
  kiwisolver=1.0.1 \
  matplotlib=3.1.3 \
  matplotlib-base=3.1.3 \
  mpc=1.1.0 \
  mpfr=4.0.2 \
  ncurses=6.1 \
  numpy=1.18.1 \
  openssl=1.1.1f \
  pandas=1.0.3 \
  pango=1.42.4 \
  pcre=8.44 \
  perl=5.26.2 \
  perl-apache-test=1.40 \
  perl-app-cpanminus=1.7044 \
  perl-archive-zip=1.64 \
  perl-autoloader=5.74 \
  perl-base=2.23 \
  perl-carp=1.38 \
  perl-class-load=0.25 \
  perl-class-load-xs=0.10 \
  perl-class-method-modifiers=2.12 \
  perl-class-std=0.013 \
  perl-clone=0.42 \
  perl-compress-raw-zlib=2.087 \
  perl-constant=1.33 \
  perl-data-dumper=2.173 \
  perl-data-optlist=0.110 \
  perl-devel-globaldestruction=0.14 \
  perl-devel-overloadinfo=0.005 \
  perl-devel-stacktrace=2.04 \
  perl-dist-checkconflicts=0.11 \
  perl-dynaloader=1.25 \
  perl-eval-closure=0.14 \
  perl-excel-writer-xlsx=1.00 \
  perl-exporter=5.72 \
  perl-extutils-makemaker=7.36 \
  perl-file-find=1.27 \
  perl-file-path=2.16 \
  perl-file-spec=3.48_01 \
  perl-file-temp=0.2304 \
  perl-getopt-long=2.50 \
  perl-http-date=6.02 \
  perl-io-tty=1.12 \
  perl-ipc-run=20180523.0 \
  perl-list-util=1.38 \
  perl-math-bigint=1.999816 \
  perl-math-complex=1.59 \
  perl-mime-base64=3.15 \
  perl-module-implementation=0.09 \
  perl-module-runtime=0.016 \
  perl-module-runtime-conflicts=0.003 \
  perl-moo=2.003004 \
  perl-moose=2.2011 \
  perl-mro-compat=0.13 \
  perl-package-deprecationmanager=0.17 \
  perl-package-stash=0.38 \
  perl-package-stash-xs=0.28 \
  perl-parallel-forkmanager=2.02 \
  perl-params-util=1.07 \
  perl-parent=0.236 \
  perl-pathtools=3.75 \
  perl-pod-escapes=1.07 \
  perl-pod-usage=1.69 \
  perl-regexp-common=2017060201 \
  perl-role-tiny=2.000008 \
  perl-scalar-list-utils=1.52 \
  perl-statistics-distributions=1.02 \
  perl-statistics-r=0.34 \
  perl-storable=3.15 \
  perl-sub-exporter=0.987 \
  perl-sub-exporter-progressive=0.001013 \
  perl-sub-identify=0.14 \
  perl-sub-install=0.928 \
  perl-sub-name=0.21 \
  perl-sub-quote=2.006003 \
  perl-test=1.26 \
  perl-test-harness=3.42 \
  perl-text-balanced=2.03 \
  perl-text-wrap=2013.0523 \
  perl-time-local=1.28 \
  perl-try-tiny=0.30 \
  perl-xml-parser=2.44 \
  perl-xsloader=0.24 \
  pip=20.0.2 \
  pixman=0.38.0 \
  pyparsing=2.4.6 \
  python=3.8.2 \
  python-dateutil=2.8.1 \
  python_abi=3.8 \
  scipy=1.4 \
  h5py=2.10 \
  pytorch torchvision $CUDA \
  pytz=2019.3 \
  r-base=3.5.1 \
  r-fnn=1.1.3 \
  r-mass=7.3_51.5 \
  readline=8.0 \
  sphinx=3.1 \
  recommonmark=0.6 
  setuptools=46.1.3 \
  six=1.14.0 \
  sqlite=3.30.1 \
  tapi=1000.10.8 \
  tk=8.6.10 \
  tktable=2.10 \
  tornado=6.0.4 \
  wheel=0.34.2 \
  xz=5.2.4 \
  zlib=1.2.11 \
  zstd=1.4.4

if [ -n "$debug" ] ; then 
    $PWD/conda/bin/perl -e 'use Config; print "using CC=$Config{cc}\n"'
fi

if [ $? -ne 0 ] ; then
    echo "Error in conda installation; aborting."
    exit 254
fi

echo "installing pip modules"
$PWD/conda/bin/pip install lxml sqlalchemy pynumpress pyteomics

if [ $? -ne 0 ] ; then 
    echo "Error in pip install; aborting."
    exit 253
fi

if [ -n "$debug" ] ; then 
    $PWD/conda/bin/perl -e 'use Config; print "using CC=$Config{cc}\n"'
fi

echo "installing cpan modules"
$PWD/conda/bin/cpanm HTTP::Message~"<= 6.20" File::Copy File::Basename Scalar::Util LWP::UserAgent Set::Partition Sys::Hostname Spreadsheet::XLSX

if [ $? -ne 0 ] ; then
    echo "Error in CPAN module installation; aborting."
    exit 252
fi

echo "configuring JUMP"
JUMP_CONFIG_PATH=$PWD/etc/cfg.bin $PWD/conda/bin/perl Makefile.PL  "PERL_BIN=$PWD/conda/bin" "$@"
if [ $? -ne 0 ] ; then
    echo "Error in JUMP configuration; aborting."
    exit 251
fi
echo "building documentation"
. $(conda env list | grep -E '^base' | awk '{print $2;}')/etc/profile.d/conda.sh
conda activate $PWD/conda
for f in {build,static,tempates} ; do if [ ! -e docs/_${f} ] ; then mkdir docs/_${f} ; fi ; done
if (cd docs && PATH=$PWD/conda/bin/python:$PATH make html) ; then if [ -e manual.html ] ; then rm manual.html ; fi ; ln -s docs/_build/html/index.html manual.html ; fi
show_success_message

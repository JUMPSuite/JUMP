use 5.6.0;
use strict;

my %options;
foreach my $arg (@ARGV) {
    if( $arg =~ /.*=.*/ ) {
	my @kv = split('=',$arg);
	my $k = shift(@kv);
	$options{$k} = join( '=', @kv );
    }
}

my $PREFIX='/usr/local/jump-g';
if(defined($options{'PREFIX'})) {
    $PREFIX=$options{'PREFIX'};
}

print "Writing Makefile for JUMPg\n";

open( my $makefile, '>Makefile' );
print $makefile <<EOF
all: 

install: 
\tif [ ! -d $PREFIX ] ; then mkdir $PREFIX ; fi
\tif [ ! -d $PREFIX/programs ] ; then mkdir $PREFIX/programs ; fi
\tfor i in \$\$(find programs -iname '*' -print) ; do if [ -d \$\$i -a ! -d $PREFIX/\$\$i ] ; then mkdir $PREFIX/\$\$i ; else if [ ! -d \$\$i ] && ([ -e $PREFIX/\$\$i ] && [ -n "\$\$(cmp \$\$i $PREFIX/\$\$i)" ]) || [ ! -e $PREFIX/\$\$i ] ; then cp -f \$\$i $PREFIX/\$\$i ; fi ; fi ; done
EOF




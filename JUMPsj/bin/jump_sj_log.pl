use Scalar::Util;
use File::Temp;

my $buf = '';
my $outf;
my @targets;

while(<STDIN>) {
    if( scalar(@targets) == 0 ) {
	$buf .= $_;
    }
    elsif (defined(Scalar::Util::openhandle($outf))) {
	print $outf $_;
    }
    else {
	open( $outf, ">$targets[0]/jump_sj.log" );
	print $outf $buf;
	print $outf $_;
	$buf = '';
    }

    if( $_ =~ m/Using: (\S+)/ ) {
	push( @targets, $1 );
    }
}

my $fname = shift( @targets );
foreach my $t (@targets) {
    link( "$fname/jump_sj.log", "$t/jump_sj.log" );
}

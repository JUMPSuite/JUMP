use File::Spec;
use Cwd;
use strict;
use Config;

sub configure_jump_driver {
    my %options = %{shift(@_)};
    my $fh;
    open( $fh, '<'.File::Spec->join('jump-driver','bin','jump.in' ) );
    my @lines;
    while(<$fh>) {
	push( @lines, $_ );
    }
    close( $fh );

    my $libpath = $Config{'installsitelib'};
    $libpath =~ s/$Config{'installprefix'}/$options{'prefix'}/g;
    
    my $header = shift(@lines);
    unshift( @lines, $header,
	     "export JUMP_READW_EXE=$options{'READW_EXE'}\n",
	     "export WINE_DIR=$options{'WINE_ROOT'}\n",
	     "PATH=".File::Spec->join($options{'PREFIX'},'bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'jump-database/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'jump-filter/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'jump-localization/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'jump-quantification/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'jump-search/bin').":\$PATH\n", 

	     "export JUMP_D_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'jump-database',$libpath))."\n", 
	     "export JUMP_F_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'jump-filter',$libpath))."\n", 
	     "export JUMP_L_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'jump-localization',$libpath))."\n", 
	     "export JUMP_Q_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'jump-quantification',$libpath))."\n", 
	     "export JUMP_SJ_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'jump-search',$libpath))."\n"
	);
    open( $fh, '>'.File::Spec->join('jump-driver','bin','jump' ) );
    print $fh @lines;
    close( $fh );
}

my %options;
foreach my $arg (@ARGV) {
    if( $arg =~ /.*=.*/ ) {
	my @kv = split('=',$arg);
	$options{$kv[0]} = $kv[1];
    }
}

unless(defined($options{'PREFIX'})) {
    $options{'PREFIX'} = $Config{'installprefix'};
}

unless(defined($options{'READW_EXE'})) {
    $options{'READW_EXE'} = 'ReAdW.exe';
}

unless(defined($options{'WINE_ROOT'})) {
    warn "without specifying WINE_ROOT, you will not be able to process RAW files.";
}

open( my $makefile, '>Makefile' );
my $optstr = join( ' ', @ARGV );
print $makefile <<EOF;
all: Makefile
\t(for i in jump-database jump-driver jump-filter jump-localization jump-quantification jump-search ; do cd \$\${i} && make && cd - ; done ; )

install: Makefile
\t(for i in jump-database jump-driver jump-filter jump-localization jump-quantification jump-search ; do cd \$\${i} && make install && cd - ; done ; )

test: Makefile
\t(for i in jump-database jump-driver jump-filter jump-localization jump-quantification jump-search ; do cd \$\${i} && make test && cd - ; done)

Makefile: Makefile.PL
\t\@echo 'Your Makefile is out-of-date with respect to Makefile.PL...reconfiguring'
\tmv Makefile Makefile.old
\t$Config{'perlpath'} Makefile.PL $optstr
\t\@echo 'Your makefiles have been rebuilt. Please re-execute your make command.'
\tfalse
EOF

my @subdirs = ('JUMPd','JUMPf', 'JUMPl', 
	       'JUMPq', 'JUMPsj');

my $rootdir = getcwd();
foreach my $dir (@subdirs) {
    chdir($dir);
    system('perl Makefile.PL PREFIX=' . File::Spec->join($options{'PREFIX'},$dir));
    chdir($rootdir);
}


configure_jump_driver(\%options);
chdir('JUMP');
system('perl Makefile.PL PREFIX=' . $options{'PREFIX'});
chdir($rootdir);


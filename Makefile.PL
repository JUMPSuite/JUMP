use File::Spec;
use Cwd;
use strict;
use Config;

sub get_etc_path {
    my %options = %{shift(@_)};
    return File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'etc'));
}

sub configure_jump_driver {
    my %options = %{shift(@_)};
    my $fh;
    open( $fh, '<'.File::Spec->join('JUMP','bin','jump.in' ) );
    my @lines;
    while(<$fh>) {
	push( @lines, $_ );
    }
    close( $fh );

    my $libpath = $Config{'installsitelib'};
    $libpath =~ s/$Config{'installprefix'}/$options{'prefix'}/g;
    
    my $header = shift(@lines);
    my @jump_g_lib_dirs = ('c','f','g','s');
    my $jump_g_root = File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPg'));
    my $etcPath = get_etc_path(\%options);
    unshift( @lines, 
	     "export JUMP_READW_EXE=$options{'READW_EXE'}\n",
	     "export PERL5LIB=\$PERL5LIB:".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'common',$libpath))."\n",
	     "PATH=".File::Spec->join($options{'PREFIX'},'bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPd/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPf/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPl/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPq/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPsj/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPi/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPv/bin').":\$PATH\n", 
	     "PATH=".File::Spec->join($options{'PREFIX'},'JUMPg/programs').":\$PATH\n", 
	     "export JUMP_D_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPd',$libpath))."\n", 
	     "export JUMP_F_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPf',$libpath))."\n", 
	     "export JUMP_L_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPl',$libpath))."\n", 
	     "export JUMP_Q_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPq',$libpath))."\n", 
	     "export JUMP_SJ_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPsj',$libpath))."\n",
	     "export JUMP_I_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPi',$libpath))."\n",
	     "export JUMP_V_LIB=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPv',$libpath))."\n",
	     "export JUMP_G_ROOT=".$jump_g_root."\n",
	     "export JUMP_CONFIG_PATH=".File::Spec->join($etcPath,'cfg.bin')."\n"
	);

    if(defined($options{'ENVIRONMENT'})) {
	while( my ($k,$v) = each %{$options{'ENVIRONMENT'}} ) {
	    unshift( @lines, "export $k=$v\n" )
	} 
    }
    unshift( @lines, $header );
    open( $fh, '>'.File::Spec->join('JUMP','bin','jump' ) );
    print $fh @lines;
    close( $fh );
}

my @subdirs = ('JUMPd','JUMPf', 'JUMPl', 
	       'JUMPq', 'JUMPsj', 'JUMPi',
	       'JUMPv', 'JUMPaq', 'JUMPg',
	       'common' );
my %kwargs;
my %options;
foreach my $arg (@ARGV) {
    if( $arg =~ /.*=.*/ ) {
	my @kv = split('=',$arg);
	my $k = shift(@kv);
	$kwargs{$k} = join( '=', @kv );
    }
}

unless(defined($kwargs{'PREFIX'})) {
    $options{'PREFIX'} = $Config{'installprefix'};
}
else {
    $options{'PREFIX'} = $kwargs{'PREFIX'};
    delete $kwargs{'PREFIX'};
}

unless(defined($kwargs{'READW_EXE'})) {
    $options{'READW_EXE'} = 'ReAdW.exe';
}
else {
    $options{'READW_EXE'} = $kwargs{'READW_EXE'};
    delete $kwargs{'READW_EXE'};
}

$options{'ENVIRONMENT'} = \%kwargs;

my $etcPath = get_etc_path(\%options);
my $cfgFile = File::Spec->join($etcPath,'cfg.bin');
my $jumpPath = File::Spec->join($options{'PREFIX'},'bin','jump');
open( my $genSiteConfigFile, '>gen-site-config.sh' );
print $genSiteConfigFile <<EOF;
if [ ! -d $etcPath ] ; then mkdir $etcPath ; fi
perl -e "use Storable; my %h; store(\\%h,'$cfgFile');"
$jumpPath -config install good
EOF

my $subdirstr = join( ' ', @subdirs ) . ' JUMP';
open( my $makefile, '>Makefile' );
my $optstr = join( ' ', @ARGV );
print $makefile <<EOF;
all: Makefile
\t(for i in $subdirstr ; do cd \$\${i} && make && cd - ; done ; )

install: Makefile
\t(for i in $subdirstr ; do cd \$\${i} && make install && cd - ; done ; sh gen-site-config.sh )

test: Makefile
\t(for i in $subdirstr ; do cd \$\${i} && make test && cd - ; done)

Makefile: Makefile.PL
\t\@echo 'Your Makefile is out-of-date with respect to Makefile.PL...reconfiguring'
\tmv Makefile Makefile.old
\t$Config{'perlpath'} Makefile.PL $optstr
\t\@echo 'Your makefiles have been rebuilt. Please re-execute your make command.'
\tfalse
EOF

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


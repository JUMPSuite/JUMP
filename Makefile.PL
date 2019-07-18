use File::Spec;
use lib File::Spec->join( 'common', 'lib' );
use Cwd;
use strict;
use Config;
use Spiders::Config;
use Storable;

sub get_etc_path {
    my %options = %{shift(@_)};
    return File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'etc'));
}

sub configure_jump_driver {
    my %options = %{shift(@_)};
    my @lines;
    my $libpath;
    my $extra_jump_dir = '';
    unless( $options{'INPLACE'} == 1 ) {
	$libpath = $Config{'installsitelib'};
	$libpath =~ s/$Config{'installprefix'}/$options{'prefix'}/g;
    }
    else {
	$libpath = 'lib';
	$extra_jump_dir = 'JUMP';
    }
    
    my $header = shift(@lines);
    my @jump_g_lib_dirs = ('c','f','g','s');
    my $jump_g_root = File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'JUMPg'));
    my $etcPath = get_etc_path(\%options);
    unshift( @lines, 
	     "export JUMP_READW_EXE=\"$options{'RAWFILE_CONVERTER'}\"\n",
	     "export PERL5LIB=\$PERL5LIB:".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},'common',$libpath))."\n",
	     "PATH=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},$extra_jump_dir,'bin')).":\$PATH\n", 
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

    open( my $fh, '>'.File::Spec->join('JUMP','bin','config.sh' ) );
    print $fh @lines;
    close( $fh );
    my $perm = (stat $fh)[2] & 07777;
    chmod( 700 | $perm, File::Spec->join('JUMP','bin','config.sh' ) );

    open( $fh, '>'.File::Spec->join('JUMP','bin','jump') );
    print $fh "#!/bin/bash\n"."PATH=".File::Spec->canonpath(File::Spec->join($options{'PREFIX'},$extra_jump_dir,'bin')).":\$PATH\n";
    open( my $template, '<'.File::Spec->join('JUMP','bin','jump.in' ) );
    while(<$template>) {
	print $fh $_;
    }
    close( $fh );
    chmod( 700 | $perm, File::Spec->join('JUMP','bin','jump' ) );
    close( $template );
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

my @defaultParamKeys = ('PREFIX',
			'RAWFILE_CONVERTER',
			'DEFAULT_DB_DIR',
			'BOOTSTRAP_DBS',
			'USE_BATCH_SYSTEM',
			'NORMAL_QUEUE',
			'LARGE_MEM_QUEUE',
			'ALLOCATION_PROJECT',
			'MAX_TIME',
			'MAX_NORMAL_MEMORY',
			'MAX_LARGE_MEMORY',
			'EXTRA_READWRITE_JOB_FLAGS',
			'EXTRA_READONLY_JOB_FLAGS',
			'INPLACE');
my @defaultParamValues = ($Config{'installprefix'},
			  'ReAdW.exe',
			  '/etc/jump',
#			  'UP000005640_9606,UP000000589_10090',
			  '',
			  'no-batch-system',
			  '',
			  '',
			  '',
			  'inf',
			  'inf',
			  'inf',
			  '',
			  '',
			  '1');
my @defaultParamHelp = ('installation prefix',
			'command for Rawfile to mzXML converter',
			'location to install default databases',
			'bootstrap databases for all FASTA files in this colon-separated list',
			'configure JUMP to use this batch system for job dispatch (e.g. PBS,LSF,slurm)',
			'queue for jobs that have no special resource needs (only if USE_BATCH_SYSTEM != no-batch)',
			'queue for jobs that require large memory (only if USE_BATCH_SYSTEM != no-batch)',
			'batch system project id to use for allocation charges (only if USE_BATCH_SYSTEM != no-batch)',
			'maximum execution time for JUMP jobs (only if USE_BATCH_SYSTEM != no-batch)',
			'maximum memory to request for jobs in the NORMAL_QUEUE (only if USE_BATCH_SYSTEM != no-batch)',
			'maximum memory to request for jobs in the LARGE_MEM_QUEUE (only if USE_BATCH_SYSTEM != no-batch)',
			'extra flags for jobs that perform intense IO',
			'extra flags for jobs that do not perform intense IO',
			'perform an in-place configuration');
my %defaultParams;
my %defaultParamHelp;
@defaultParams{@defaultParamKeys} = @defaultParamValues;
@defaultParamHelp{@defaultParamKeys} = @defaultParamHelp;

if( grep( /^-h$/, @ARGV ) ) {
    print <<EOF;
usage: perl Makefile.PL <configuration key=value pairs>
\tconfiguration script for JUMP suite
    
Some influmential configuration key=value pairs:

EOF

foreach my $k (@defaultParamKeys) {
    my $v = $defaultParams{$k};
    if( $v eq '' ) {
	print "\t$k=$defaultParamHelp{$k}\n\t\tno default value\n\n";
    }
    else {
	print "\t$k=$defaultParamHelp{$k}\n\t\tdefault:$v\n\n";
    }
}
    print <<EOF;
You may pass any key=value pair to Makefile.PL.  Any key=value pair not
listed above will be added as an environment variable for JUMP to use.
EOF

    exit 1;
}

while( my ($k,$v) = each %defaultParams ) {
    unless(defined($kwargs{$k})) {
	$options{$k} = $defaultParams{$k};
    }
    else {
	$options{$k} = $kwargs{$k};
	delete $kwargs{$k};
    }
}

my $rootdir = $ENV{"PWD"};
if( $options{'INPLACE'} == 1 ) {
    $options{'PREFIX'} = $rootdir;
}

if( $options{'USE_BATCH_SYSTEM'} ne 'no-batch-system' ) {
    my $batch_complete = 1;
    foreach my $k ('NORMAL_QUEUE','LARGE_MEM_QUEUE','ALLOCATION_PROJECT',
		   'MAX_TIME','MAX_NORMAL_MEMORY','MAX_LARGE_MEMORY') {
	if($options{$k} eq '') {
	    print "You must define $k since you have selected USE_BATCH_SYSTEM=$options{'USE_BATCH_SYSTEM'}\n";
	    $batch_complete = 0;
	}
    }
    if( ! $batch_complete ) {
	exit -1;
    }
}

$options{'CLUSTER'} = ($options{'USE_BATCH_SYSTEM'} ne 'no-batch-system');
my $etcPath = get_etc_path(\%options);
my $cfgFile = File::Spec->join($etcPath,'cfg.bin');
if( ! -d $etcPath ) {
    mkdir( $etcPath )
}

{my %h; store(\%h,"$cfgFile");}
$ENV{'JUMP_CONFIG_PATH'} = $cfgFile;
my $config = new Spiders::Config();
$config->put( "use_batch_system", $options{'USE_BATCH_SYSTEM'} );
$config->put( "normal_queue", $options{'NORMAL_QUEUE'} );
$config->put( "large_mem_queue", $options{'LARGE_MEM_QUEUE'} );
$config->put( "allocation_project", $options{'ALLOCATION_PROJECT'} );
$config->put( "max_time", $options{'MAX_TIME'} );
$config->put( "max_normal_memory", $options{'MAX_NORMAL_MEMORY'} );
$config->put( "max_large_memory", $options{'MAX_LARGE_MEMORY'} );
$config->put( "extra_readwrite_job_flags", $options{'EXTRA_READWRITE_JOB_FLAGS'} );
$config->put( "extra_readonly_job_flags", $options{'EXTRA_READONLY_JOB_FLAGS'} );
$ENV{'JUMP_CONFIG_PATH'} = '';

$options{'ENVIRONMENT'} = \%kwargs;

my $subdirstr = join( ' ', @subdirs ) . ' JUMP';


if($options{'INPLACE'} eq '1') {
    configure_jump_driver(\%options);
}
else {
    foreach my $dir (@subdirs) {
	chdir($dir);
	system('perl Makefile.PL PREFIX=' . File::Spec->join($options{'PREFIX'},$dir));
	chdir($rootdir);
    }

    configure_jump_driver(\%options);
    chdir('JUMP');
    system('perl Makefile.PL PREFIX=' . $options{'PREFIX'});
    chdir($rootdir);

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
}


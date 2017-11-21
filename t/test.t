use File::Copy;
use File::Path;
use File::Basename;

use Test::Simple tests => 4;

my $param_file = "example/jump.params";
my $abs_param_file = File::Spec->rel2abs( $param_file );
my $mz_file = "example/HH_tmt10_human_jump_truncated2048.mzXML";
my $abs_mz_file = File::Spec->rel2abs( $mz_file );
my $mz_root = $mz_file;
$mz_root =~ s/.mzXML//;

# relative everything
system("perl bin/jump_sj.pl -p $param_file $mz_file &> /dev/null");
ok( !$?, "relative params, relative mzXML" );
move( File::Spec->join( $mz_root, basename($mz_file) ), $mz_file );
File::Path->remove_tree( $mz_root );

# absolute MZ file
system("perl bin/jump_sj.pl -p $param_file $abs_mz_file &> /dev/null");
ok( !$?, "relative params, absolute mzXML" );
move( File::Spec->join( $mz_root, basename($mz_file) ), $mz_file );
File::Path->remove_tree( $mz_root );

# absolute param file
system("perl bin/jump_sj.pl -p $abs_param_file $mz_file &> /dev/null");
ok( !$?, "absolute params, relative mzXML" );
move( File::Spec->join( $mz_root, basename($mz_file) ), $mz_file );
File::Path->remove_tree( $mz_root );

# absolute everything
system("perl bin/jump_sj.pl -p $abs_param_file $abs_mz_file &> /dev/null");
ok( !$?, "absolute params, absolute mzXML" );
move( File::Spec->join( $mz_root, basename($mz_file) ), $mz_file );
File::Path->remove_tree( $mz_root );


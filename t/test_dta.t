use Spiders::Dta;
use Spiders::Dtas;
use Spiders::FsDtasBackend;
use Spiders::IdxDtasBackend;

$fs_dtas = Spiders::Dtas->new(Spiders::FsDtasBackend->new("HH_tmt10_human_jump_truncated2048/HH_tmt10_human_jump_truncated2048.72","read"));
$idx_dtas = Spiders::Dtas->new(Spiders::IdxDtasBackend->new("HH_tmt10_human_jump_truncated2048/HH_tmt10_human_jump_truncated2048.75","read"));

for $k (@{$fs_dtas->list_dta()}) {
    $fsdta = $fs_dtas->get_dta($k);
    $idxdta = $idx_dtas->get_dta($k);
    $fshash = $fsdta->get_mz_int_hash();
    $idxhash = $idxdta->get_mz_int_hash();
    @fskeys = keys(%$fshash);
    @idxkeys = keys(%$idxhash);
    @fskeys = sort {$a <=> $b} @fskeys;
    @idxkeys = sort {$a <=> $b} @idxkeys;
    $fsmz = $fsdta->get_mz_array();
    $idxmz = $idxdta->get_mz_array();
    $fsint = $fsdta->get_int_array();
    $idxint = $idxdta->get_int_array();
    for( $i = 0; $i < scalar(@$fsmz); ++$i ) {
	$mzfs = $$fsmz[$i];
	$mzidx = $$idxmz[$i];
	$intfs = $$fsint[$i];
	$intidx = $$idxint[$i];
	if( abs($mzfs - $mzidx) > 1e-5 ) {
	    print "mz error in file ",$k,"\n";
	    last;
	}
	if( abs($intfs - $intidx) > 1e-5 ) {
	    print "int error in file ",$k,"\n";
	    last;
	}
	if( abs($fskeys[$i] - $mzfs) > 1e-5 ||
	    abs($idxkeys[$i] - $mzidx) > 1e-5 ||
	    abs($fshash->{$fskeys[$i]} - $idxhash->{$idxkeys[$i]}) > 1e-5 ) {
	    print "hash error in file $k at index $i\n";
	    last;
	}
    }
}

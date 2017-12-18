#!/usr/bin/perl
#
package Spiders::WilcoxonRankSum;

use strict;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA     = qw(Exporter);
@EXPORT      = qw(_check_dataset  _compute_N_MaxSum  _reset_dependant_datastructures  _rank_sum_for  _set_smaller_rank_for  _NormalZ  set_dataset1  set_dataset2  load_data  compute_ranks  compute_rank_array  rank_sum_for  get_smaller_rank_sum  smaller_rank_sums_count  rank_sums_other_than_expected_counts  probability_exact  bcomb  probability_normal_approx  probability_status  chisqrdistr  udistr  tdistr  fdistr  uprob  chisqrprob  tprob  fprob  _subfprob  _subchisqrprob  _subu  _subuprob  _subt  _subtprob  _subf  _subf2  _subchisqr  log10  precision  precision_string);

use List::Util qw(sum);
use Set::Partition;
use Math::BigInt;


use Class::Std;

sub new
{
    my ($class,%arg) = @_;
    my $self = {
    };
    bless $self, $class;
     return $self;
}

use constant SIGNIFICANT => 5; 

use constant PI => 3.14159265;

{
    my %EXACT_UPTO : ATTR( :init_arg<exact_upto> :default<20> );
    my %dataset1_of : ATTR( :get<dataset1> );   
    my %dataset2_of : ATTR( :get<dataset2> );   
    my %n1_of : ATTR( :get<n1> );               
    my %n2_of : ATTR( :get<n2> );               
    my %N_of : ATTR( :get<N>  );     
    my %MaxSum_of : ATTR( :get<max_rank_sum> );     
    my %ranks_of : ATTR( :get<ranks>    :set<ranks>  );  
    my %rank_array_of : ATTR( :get<rank_array> );        
    my %rankSum1_of : ATTR( :get<rankSum_dataset1> );    
    my %expected_rank_sum_1_of : ATTR( :get<expected_rank_sum_dataset1>)
      ;     
    my %expected_rank_sum_2_of : ATTR( :get<expected_rank_sum_dataset2>)
      ;     
    my %rankSum2_of : ATTR( :get<rankSum_dataset2> );    
    my %smaller_rank_sum_of : ATTR;
    my %smaller_ranks_count_of : ATTR;
    my %expected_rank_count_for_smaller_ranks_count_of :
      ATTR( :get<expected_rank_count_for_smaller_ranks_count>);
    my %smaller_rank_sums_count_of :
      ATTR;     
                
    my %rank_sums_other_than_expected_count_of :
      ATTR;     
                
    my %probability_of :
      ATTR;     
    my %probability_normal_approx_of : ATTR;

     

    sub _check_dataset {
        my ($dataset_ref) = @_;

        die "Need array ref to dataset\n" unless ($dataset_ref);

        die "Datasets must be passed as array references\n" unless ( ref($dataset_ref) eq 'ARRAY' );

        my @dataset = grep { $_ > 0 } @{$dataset_ref};
        die "dataset has no element greater 0\n" unless (@dataset);

        return \@dataset;

    }

    sub _compute_N_MaxSum {
        my ($id) = @_;

        my $N;
        unless ( $N_of{$id} ) {
            $N = $n1_of{$id} + $n2_of{$id};
            $N_of{$id} = $N;
        }

        unless ( $MaxSum_of{$id} ) {
            $MaxSum_of{$id} = $N * ( $N + 1 ) / 2;
        }

        unless ( $expected_rank_sum_1_of{$id} ) {
            $expected_rank_sum_1_of{$id} = $n1_of{$id} * $N / 2;
        }

        unless ( $expected_rank_sum_2_of{$id} ) {
            $expected_rank_sum_2_of{$id} = $n2_of{$id} * $N / 2;
        }

        return;
    }

    sub _reset_dependant_datastructures {
        my ($id) = @_;

        delete $ranks_of{$id};
        delete $rank_array_of{$id};
        delete $rankSum1_of{$id};
        delete $rankSum2_of{$id};
        delete $N_of{$id};
        delete $MaxSum_of{$id};
        delete $smaller_rank_sum_of{$id};
        delete $smaller_rank_sums_count_of{$id};
        delete $probability_of{$id};
        delete $probability_normal_approx_of{$id};
        delete $expected_rank_sum_1_of{$id};
        delete $expected_rank_sum_2_of{$id};

        return;
    }

    sub _rank_sum_for {
        my ( $self, $dataset ) = @_;

        my $id = ident $self;

        my @rank_array;

        if ( $rank_array_of{$id} and @{ $rank_array_of{$id} } ) {
            @rank_array = @{ $rank_array_of{$id} };
        }
        else {
            @rank_array = $self->compute_rank_array();
        }

        return sum map { $_->[0] } grep { $_->[1] eq $dataset } @rank_array;
    }

    sub _set_smaller_rank_for {
        my ( $id, $rank_sum_1, $rank_sum_2 ) = @_;
        if ( $rank_sum_1 <= $rank_sum_2 ) {
            $smaller_rank_sum_of{$id}    = $rank_sum_1;
            $smaller_ranks_count_of{$id} = $n1_of{$id};
            $expected_rank_count_for_smaller_ranks_count_of{$id} =
              $expected_rank_sum_1_of{$id};
        }
        else {
            $smaller_rank_sum_of{$id}    = $rank_sum_2;
            $smaller_ranks_count_of{$id} = $n2_of{$id};
            $expected_rank_count_for_smaller_ranks_count_of{$id} =
              $expected_rank_sum_2_of{$id};
        }
        return;
    }

    sub _NormalZ {     
        my ($x) = @_;
         
         
         
         
         
         
        my @b =
          ( 0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429 );
        my $p = 0.2316419;
        my $t = 1 / ( 1 + $p * $x );

         
        my $fact = $t;
        my $Sum;

         
        foreach my $bi (@b) {
            $Sum += $bi * $fact;
            $fact *= $t;
        }

         
        $p =
          2 * $Sum *
          exp( -$x * $x / 2.0 ) /
          ( sqrt( 2 * 3.14159265358979323846 ) );
         
        return $p;
    }

     

    sub set_dataset1 {
        my ( $self, $dataset1_ref ) = @_;

        $dataset1_ref = _check_dataset($dataset1_ref);

        my $id = ident $self;
        $dataset1_of{$id} = $dataset1_ref;
        $n1_of{$id}       = scalar( @{$dataset1_ref} );

        _reset_dependant_datastructures($id);

        return;
    }

    sub set_dataset2 {
        my ( $self, $dataset2_ref ) = @_;

        $dataset2_ref = _check_dataset($dataset2_ref);

        my $id = ident $self;
        $dataset1_of{$id} = $dataset2_ref;
        $n2_of{$id}       = scalar( @{$dataset2_ref} );

        _reset_dependant_datastructures($id);

        return;
    }

    sub load_data {
        my ( $self, $dataset1_ref, $dataset2_ref ) = @_;

        $dataset1_ref = _check_dataset($dataset1_ref);
        $dataset2_ref = _check_dataset($dataset2_ref);

        my $id = ident $self;

        $dataset1_of{$id} = $dataset1_ref;
        $dataset2_of{$id} = $dataset2_ref;
        $n1_of{$id}       = scalar( @{$dataset1_ref} );
        $n2_of{$id}       = scalar( @{$dataset2_ref} );

        _reset_dependant_datastructures($id);

        _compute_N_MaxSum($id);

        return;
    }

    sub compute_ranks {
        my ($self) = @_;
        my $id = ident $self;

        die "Please set/load datasets before computing ranks\n"
          unless ( $dataset1_of{$id} and $dataset2_of{$id} );

        my @dataset1 = @{ $dataset1_of{$id} };
        my @dataset2 = @{ $dataset2_of{$id} };

 
        _compute_N_MaxSum($id);

        my %ranks;

        foreach my $el (@dataset1) {
            $ranks{$el}->{in_dataset}->{ds1}++;
        }
        foreach my $el (@dataset2) {
            $ranks{$el}->{in_dataset}->{ds2}++;
        }

        my $rank = 0;
        foreach my $value ( sort { $a <=> $b } keys %ranks ) {

            my $tied_ranks;

            foreach my $ds ( keys %{ $ranks{$value}->{in_dataset} } ) {
                $tied_ranks += $ranks{$value}->{in_dataset}->{$ds};
            }

 

            my $rs;
            for my $r ( $rank + 1 .. $rank + $tied_ranks ) {
                $rs += $r;
            }
            $ranks{$value}->{rank} = $rs / $tied_ranks;
            $ranks{$value}->{tied} = $tied_ranks;

            $rank += $tied_ranks;
        }

        $ranks_of{$id} = \%ranks;

        return $ranks_of{$id};
    }

    sub compute_rank_array {
        my ($self) = @_;
        my $id = ident $self;

        my @rank_array;
        if ( $rank_array_of{$id} and @{ $rank_array_of{$id} } ) {
            @rank_array = @{ $rank_array_of{$id} };
        }
        else {

            my %ranks;

            if ( $ranks_of{$id} and %{ $ranks_of{$id} } ) {
                %ranks = %{ $ranks_of{$id} };
            }
            else {
                %ranks = %{ $self->compute_ranks() };
            }

            foreach my $value ( sort { $a <=> $b } keys %ranks ) {
                foreach my $ds ( keys %{ $ranks{$value}->{in_dataset} } ) {
                    for ( 1 .. $ranks{$value}->{in_dataset}->{$ds} ) {
                        push( @rank_array, [ $ranks{$value}->{rank}, $ds ] );
                    }
                }
            }

            $rank_array_of{$id} = \@rank_array;

        }

        return (@rank_array);
    }

    sub rank_sum_for {
        my ( $self, $for_dataset ) = @_;

        my $id = ident $self;

        my $rankSum;
        if ( $for_dataset =~ m{1} ) {
            if ( $rankSum1_of{$id} ) {
                return $rankSum1_of{$id};
            }
            else {
                $rankSum1_of{$id} = $self->_rank_sum_for('ds1');
                return $rankSum1_of{$id};
            }
        }
        elsif ( $for_dataset =~ m{2} ) {
            if ( $rankSum2_of{$id} ) {
                return $rankSum2_of{$id};
            }
            else {
                $rankSum2_of{$id} = $self->_rank_sum_for('ds2');
                return $rankSum2_of{$id};
            }
        }
        else {
            die "Argument must match `1' or `2' (meaning dataset 1 or 2)\n";
        }

        return;

    }

    sub get_smaller_rank_sum {
        my ($self) = @_;

        my $id = ident $self;

        if ( $smaller_rank_sum_of{$id} and $smaller_ranks_count_of{$id} ) {

            return (
                    ( $smaller_rank_sum_of{$id}, $smaller_ranks_count_of{$id} )
            );
        }

        my $rank_sum_1 = $rankSum1_of{$id};
        my $rank_sum_2 = $rankSum2_of{$id};

        if ( not($rank_sum_1) and not($rank_sum_2) ) {
            $rank_sum_1 = $self->rank_sum_for('ds1');
        }

        if ( $rank_sum_1 and $rank_sum_2 ) {

            _set_smaller_rank_for( $id, $rank_sum_1, $rank_sum_2 );

        }
        elsif ($rank_sum_1) {
            $rank_sum_2 = $MaxSum_of{$id} - $rank_sum_1;
            $rankSum2_of{$id} = $rank_sum_2;

            _set_smaller_rank_for( $id, $rank_sum_1, $rank_sum_2 );

        }
        elsif ($rank_sum_2) {
            $rank_sum_1 = $MaxSum_of{$id} - $rank_sum_2;
            $rankSum1_of{$id} = $rank_sum_1;

            _set_smaller_rank_for( $id, $rank_sum_1, $rank_sum_2 );

        }

        return (
                ( $smaller_rank_sum_of{$id}, $smaller_ranks_count_of{$id} )
        );

        return $smaller_rank_sum_of{$id};
    }

    sub smaller_rank_sums_count {
        my ($self) = @_;
        my $id = ident $self;

        if ( $smaller_rank_sums_count_of{$id} ) {
            return $smaller_rank_sums_count_of{$id};
        }

        my ( $W, $nA ) = $self->get_smaller_rank_sum();
        my $N      = $N_of{$id};
        my $nB     = $N - $nA;
        my $MaxSum = $MaxSum_of{$id};

        my @ranks = map { $_->[0] } $self->compute_rank_array();

         
        if ( $W > $MaxSum ) {
            die "Rank sum bound $W is bigger than the maximum possible rank sum $MaxSum\n";
        }
        if ( $N != scalar(@ranks) ) {
            die "Sum of $nA and $nB must be equal to number of ranks: "
              . scalar(@ranks) . "\n";
        }

         
        my $s = Set::Partition->new(
            list      => \@ranks,
            partition => [ $nA, $nB ],
        );

        my $count_less_W = 0;

        while ( my $p = $s->next() ) {
            my @pA   = @{ $p->[0] };
            my $sumA = sum @pA;
            if ( $sumA <= $W ) {
                $count_less_W++;
            }
        }

        return $count_less_W;

    }

    sub rank_sums_other_than_expected_counts {
        my ($self) = @_;
        my $id = ident $self;

        if ( $rank_sums_other_than_expected_count_of{$id} ) {
            return $rank_sums_other_than_expected_count_of{$id};
        }

        my ( $W, $nA ) = $self->get_smaller_rank_sum();
        my $W_exp = $self->get_expected_rank_count_for_smaller_ranks_count();

        my $N      = $N_of{$id};
        my $nB     = $N - $nA;
        my $MaxSum = $MaxSum_of{$id};

        my @ranks = map { $_->[0] } $self->compute_rank_array();

         
        if ( $W > $MaxSum ) {
            die
"Rank sum bound $W is bigger than the maximum possible rank sum $MaxSum\n";
        }
        if ( $N != scalar(@ranks) ) {
            die "Sum of $nA and $nB must be equal to number of ranks: "
              . scalar(@ranks) . "\n";
        }

         
        my $s = Set::Partition->new(
            list      => \@ranks,
            partition => [ $nA, $nB ],
        );

        my $count_other_W = 0;

        if ( $W >= $W_exp ) {

            while ( my $p = $s->next() ) {
                my @pA   = @{ $p->[0] };
                my $sumA = sum @pA;
                if ( $sumA >= $W ) {
                    $count_other_W++;
                }
            }

        }
        else {

            while ( my $p = $s->next() ) {
                my @pA   = @{ $p->[0] };
                my $sumA = sum @pA;
                if ( $sumA <= $W ) {
                    $count_other_W++;
                }
            }

        }

        return $count_other_W;

    }

    sub probability : NUMERIFY {
        my ($self) = @_;
        my $id = ident $self;

        if ( $probability_of{$id} ) {
            return $probability_of{$id};
        }

        my ( $W, $nA ) = $self->get_smaller_rank_sum();
        my $N = $N_of{$id};

        my $p;
     
      
       
        
            $p = $self->probability_normal_approx();
         

        $probability_of{$id} = $p;

        return $probability_of{$id};
    }

    sub probability_exact {
        my ($self) = @_;
        my $id = ident $self;

        my ( $W, $nA ) = $self->get_smaller_rank_sum();
        my $N = $N_of{$id};

        my $partition_count = bcomb( $N, $nA );
        my $have_smaller_rank_sums =
          $self->rank_sums_other_than_expected_counts();
       my $p =
          Math::BigFloat->new($have_smaller_rank_sums) * 2.0 /
          Math::BigFloat->new($partition_count);
        if ( $p > 1 ) { $p = 1 }

        return $p;
    }
	
	sub bcomb {
    my ($self) = shift @_;	
    my( $n, $k ) = @_;
    $n = Math::BigInt->new( $n );
    $k = Math::BigInt->new( $k );
    my $r = $n - $k;
    return $n->bfac() / ($k->bfac() * $r->bfac());
	}


    sub probability_normal_approx {
        my ($self) = @_;
        my $id = ident $self;

        my ( $W, $nA ) = $self->get_smaller_rank_sum();
        my $N          = $N_of{$id};
        my $nB         = $N - $nA;
        my $mean       = $nA * ( $N + 1 ) / 2;
        my $deviation  = sqrt( $nA * $nB * ( $N + 1 ) / 12.0 );
        my $continuity = ( ( $W - $mean ) >= 0 ) ? -0.5 : +0.5;
        my $z          = ( $W - $mean + $continuity ) / $deviation;
        @{ $probability_normal_approx_of{$id} }{ 'mean', 'std deviation', 'z' }
          = ( $mean, $deviation, $z );
 
 
 
 
 
		my $p = uprob( abs($z) );
 
		if($W<($mean) )
		{
			$p = 1 - $p;
		}
 
        return $p;

    }

    sub probability_status {
        my ($self) = (@_);
        my $id = ident $self;

        my $return_string;
        if ( $probability_of{$id} ) {
            if ( $probability_normal_approx_of{$id} ) {
                $return_string =
                  sprintf
"Probability: %10f, normal approx w. mean: %10f, std deviation: %10f, z: %10f",
                  $probability_of{$id},
                  map { $probability_normal_approx_of{$id}->{$_} }
                  ( 'mean', 'std deviation', 'z' );
            }
            else {
                $return_string = sprintf "Probability: %10f, exact",
                  $probability_of{$id};
            }
        }
        else {
            $return_string = "Probability not yet computed";
        }

        return ( "$return_string");
    }

    sub as_hash : HASHIFY {
        my ($self) = @_;
        my $id = ident $self;

        return {
            dataset_1                 => $dataset1_of{$id},
            dataset_2                 => $dataset2_of{$id},
            n1                        => $n1_of{$id},
            n2                        => $n2_of{$id},
            N                         => $N_of{$id},
            rank_array                => $rank_array_of{$id},
            rank_sum_1                => $rankSum1_of{$id},
            rank_sum_2                => $rankSum2_of{$id},
            rank_sum_1_expected       => $expected_rank_sum_1_of{$id},
            rank_sum_2_expected       => $expected_rank_sum_2_of{$id},
            probability               => $probability_of{$id},
            probability_normal_approx => $probability_normal_approx_of{$id},
        };

    }

    sub summary : STRINGIFY {
        my ($self) = (@_);
        my $id = ident $self;

        my $hash = $self->as_hash();

        my $return_string;
        if ( not( $hash->{dataset_1} ) ) {
            $return_string =
"Dataset 1 is not yet initialised, no computations could be done\n";
        }
        elsif ( not( $hash->{dataset_2} ) ) {
            $return_string =
"Dataset 2 is not yet initialised, no computations could be done\n";
        }
        else {
            my $format = <<END_FORMAT;
----------------------------------------------------------------
dataset |    n      | rank sum: observed / expected 
----------------------------------------------------------------
   1    |%7d    |           %7d      /%7d
----------------------------------------------------------------
   2    |%7d    |           %7d      /%7d
----------------------------------------------------------------
N (size of both datasets): %7d
%s
END_FORMAT
            my $prob = $self->probability_status();
            $return_string = sprintf $format, @{$hash}{
                'n1', 'rank_sum_1', 'rank_sum_1_expected', 'n2', 'rank_sum_2',
                'rank_sum_2_expected', 'N'
              },
              $prob;
            if ( $hash->{probability} >= 0.05 ) {
                $return_string .= "Not significant (at 0.05 level)\n";
            }
            else {
                $return_string .= "Significant (at 0.05 level)\n";
                $return_string .=
                  $hash->{rank_sum_1} > $hash->{rank_sum_1_expected}
                  ? "Ranks of dataset 1 are higher than expected\n"
                  : "Ranks of dataset 1 are lower than expected\n";

            }
            if ( $hash->{N} < 5 ) {
                $return_string .=
                  "Warning: sample size ($hash->{N}) too small (<5)!\n";
            }
        }

        return ( STR { "$return_string" } VOID { print $return_string } );
    }

}

sub chisqrdistr {  
		my ($n, $p) = @_;
			if ($n <= 0 || abs($n) - abs(int($n)) != 0) {
						die "Invalid n: $n\n";  
								}
				if ($p <= 0 || $p > 1) {
							die "Invalid p: $p\n"; 
								}
					return precision_string(_subchisqr($n, $p));
}

sub udistr {  
		my ($p) = (@_);
			if ($p > 1 || $p <= 0) {
						die "Invalid p: $p\n";
							}
				return precision_string(_subu($p));
}

sub tdistr {  
		my ($n, $p) = @_;
			if ($n <= 0 || abs($n) - abs(int($n)) != 0) {
						die "Invalid n: $n\n";
							}
				if ($p <= 0 || $p >= 1) {
							die "Invalid p: $p\n";
								}
					return precision_string(_subt($n, $p));
}

sub fdistr {  
		my ($n, $m, $p) = @_;
			if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
						die "Invalid n: $n\n";  
								}
				if (($m<=0) || ((abs($m)-(abs(int($m))))!=0)) {
							die "Invalid m: $m\n";  
									}
					if (($p<=0) || ($p>1)) {
								die "Invalid p: $p\n";
									}
						return precision_string(_subf($n, $m, $p));
}

sub uprob {  
		my ($x) = @_;
			return precision_string(_subuprob($x));
}

sub chisqrprob {  
		my ($n,$x) = @_;
			if (($n <= 0) || ((abs($n) - (abs(int($n)))) != 0)) {
						die "Invalid n: $n\n";  
								}
				return precision_string(_subchisqrprob($n, $x));
}

sub tprob {  
		my ($n, $x) = @_;
			if (($n <= 0) || ((abs($n) - abs(int($n))) !=0)) {
						die "Invalid n: $n\n";  
								}
				return precision_string(_subtprob($n, $x));
}

sub fprob {  
		my ($n, $m, $x) = @_;
			if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
						die "Invalid n: $n\n";  
								}
				if (($m<=0) || ((abs($m)-(abs(int($m))))!=0)) {
							die "Invalid m: $m\n";  
									} 
					return precision_string(_subfprob($n, $m, $x));
}


sub _subfprob {
	my ($n, $m, $x) = @_;
	my $p;
	if ($x<=0) 
	{
		$p=1;
	} 
	elsif ($m % 2 == 0)
	{
		my $z = $m / ($m + $n * $x);
		my $a = 1;
		for (my $i = $m - 2; $i >= 2; $i -= 2) 
		{
			$a = 1 + ($n + $i - 2) / $i * $z * $a;
		}
		$p = 1 - ((1 - $z) ** ($n / 2) * $a);
	}
	elsif ($n % 2 == 0) {
		my $z = $n * $x / ($m + $n * $x);
		my $a = 1;
		for (my $i = $n - 2; $i >= 2; $i -= 2) {
			$a = 1 + ($m + $i - 2) / $i * $z * $a;
		}
		$p = (1 - $z) ** ($m / 2) * $a;
	} 
	else 
	{
		my $y = atan2(sqrt($n * $x / $m), 1);
		my $z = sin($y) ** 2;
		my $a = ($n == 1) ? 0 : 1;
		for (my $i = $n - 2; $i >= 3; $i -= 2) {
		$a = 1 + ($m + $i - 2) / $i * $z * $a;
		} 
		my $b = PI;
		for (my $i = 2; $i <= $m - 1; $i += 2) {
			$b *= ($i - 1) / $i;
		}
		my $p1 = 2 / $b * sin($y) * cos($y) ** $m * $a;
			
		$z = cos($y) ** 2;
		$a = ($m == 1) ? 0 : 1;
		for (my $i = $m-2; $i >= 3; $i -= 2) {
			$a = 1 + ($i - 1) / $i * $z * $a;
		}
		$p = max(0, $p1 + 1 - 2 * $y / PI - 2 / PI * sin($y) * cos($y) * $a);
	}
	return $p;
}


sub _subchisqrprob {
	my ($n,$x) = @_;
	my $p;
	if ($x <= 0) {
		$p = 1;
	}
	elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)	- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} 
	elsif ($x > 400) {
		$p = 0;
	}
	else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
			$i1 = 1;
		}
		else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}
		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
}

sub _subu {
	my ($p) = @_;
	my $y = -log(4 * $p * (1 - $p));
	my $x = sqrt($y * (1.570796288 + $y * (.03706987906 + $y * (-.8364353589E-3	+ $y *(-.2250947176E-3 + $y * (.6841218299E-5 + $y * (0.5824238515E-5+ $y * (-.104527497E-5	+ $y * (.8360937017E-7+ 
	$y * (-.3231081277E-8+ $y * (.3657763036E-10+ $y *.6936233982E-12)))))))))));
	$x = -$x if ($p>.5);
	return $x;
}

sub _subuprob {
	my ($x) = @_;
	my $p = 0;  
	my $absx = abs($x);
	if ($absx < 1.9) {
	$p = (1 +
	$absx * (.049867347
	+ $absx * (.0211410061
	+ $absx * (.0032776263
	+ $absx * (.0000380036
	+ $absx * (.0000488906
	+ $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
	for (my $i = 18; $i >= 1; $i--) {
	$p = $i / ($absx + $p);
	}
	$p = exp(-.5 * $absx * $absx) 
	/ sqrt(2 * PI) / ($absx + $p);
	}
	$p = 1 - $p if ($x<0);
	return $p;
}

   
sub _subt {
	my ($n, $p) = @_;
	if ($p >= 1 || $p <= 0) {
	die "Invalid p: $p\n";
	}
	if ($p == 0.5) {
	return 0;
	} elsif ($p < 0.5) {
	return - _subt($n, 1 - $p);
	}
	my $u = _subu($p);
	my $u2 = $u ** 2;
	my $a = ($u2 + 1) / 4;
	my $b = ((5 * $u2 + 16) * $u2 + 3) / 96;
	my $c = (((3 * $u2 + 19) * $u2 + 17) * $u2 - 15) / 384;
	my $d = ((((79 * $u2 + 776) * $u2 + 1482) * $u2 - 1920) * $u2 - 945) / 92160;
	my $e = (((((27 * $u2 + 339) * $u2 + 930) * $u2 - 1782) * $u2 - 765) * $u2 + 17955) / 368640;
	my $x = $u * (1 + ($a + ($b + ($c + ($d + $e / $n) / $n) / $n) / $n) / $n);
	
	if ($n <= log10($p) ** 2 + 3) {
	my $round;
	do { 
	my $p1 = _subtprob($n, $x);
	my $n1 = $n + 1;
	my $delta = ($p1 - $p) 
	/ exp(($n1 * log($n1 / ($n + $x * $x)) 
	+ log($n/$n1/2/PI) - 1 
	+ (1/$n1 - 1/$n) / 6) / 2);
	$x += $delta;
	$round = sprintf("%.".abs(int(log10(abs $x)-4))."f",$delta);
	} while (($x) && ($round != 0));
	}
	return $x;
}

sub _subtprob {
	my ($n, $x) = @_;

	my ($a,$b);
	my $w = atan2($x / sqrt($n), 1);
	my $z = cos($w) ** 2;
	my $y = 1;
	
	for (my $i = $n-2; $i >= 2; $i -= 2) {
	$y = 1 + ($i-1) / $i * $z * $y;
	} 
	if ($n % 2 == 0) {
	$a = sin($w)/2;
	$b = .5;
	} else {
	$a = ($n == 1) ? 0 : sin($w)*cos($w)/PI;
	$b= .5 + $w/PI;
	}
	return max(0, 1 - $b - $a * $y);
}

sub _subf {
	my ($n, $m, $p) = @_;
	my $x;

	if ($p >= 1 || $p <= 0) {
	die "Invalid p: $p\n";
	}
	if ($p == 1) {
	$x = 0;
	} elsif ($m == 1) {
	$x = 1 / (_subt($n, 0.5 - $p / 2) ** 2);
	} elsif ($n == 1) {
	$x = _subt($m, $p/2) ** 2;
	} elsif ($m == 2) {
	my $u = _subchisqr($m, 1 - $p);
	my $a = $m - 2;
	$x = 1 / ($u / $m * (1 +
	(($u - $a) / 2 +
	(((4 * $u - 11 * $a) * $u + $a * (7 * $m - 10)) / 24 +
	(((2 * $u - 10 * $a) * $u + $a * (17 * $m - 26)) * $u
	- $a * $a * (9 * $m - 6)
	)/48/$n
	)/$n
	)/$n));
	} elsif ($n > $m) {
	$x = 1 / _subf2($m, $n, 1 - $p)
	} else {
	$x = _subf2($n, $m, $p)
	}
	return $x;
}

sub _subf2 {
	my ($n, $m, $p) = @_;
	my $u = _subchisqr($n, $p);
	my $n2 = $n - 2;
	my $x = $u / $n *
	(1 + 
	(($u - $n2) / 2 + 
	(((4 * $u - 11 * $n2) * $u + $n2 * (7 * $n - 10)) / 24 + 
	(((2 * $u - 10 * $n2) * $u + $n2 * (17 * $n - 26)) * $u 
	- $n2 * $n2 * (9 * $n - 6)) / 48 / $m) / $m) / $m);
	my $delta;
	do {
		my $z = exp(
		(($n+$m) * log(($n+$m) / ($n * $x + $m)) + ($n - 2) * log($x)
		+ log($n * $m / ($n+$m))- log(4 * PI)- (1/$n  + 1/$m - 1/($n+$m))/6	)/2);
		$delta = (_subfprob($n, $m, $x) - $p) / $z;
		$x += $delta;
	} while (abs($delta)>3e-4);
	
	return $x;
}

sub _subchisqr {
	my ($n, $p) = @_;
	my $x;
	
	if (($p > 1) || ($p <= 0)) {
	die "Invalid p: $p\n";
	} elsif ($p == 1){
	$x = 0;
	} elsif ($n == 1) {
	$x = _subu($p / 2) ** 2;
	} elsif ($n == 2) {
	$x = -2 * log($p);
	} else {
	my $u = _subu($p);
	my $u2 = $u * $u;
	
	$x = max(0, $n + sqrt(2 * $n) * $u 
	+ 2/3 * ($u2 - 1)
	+ $u * ($u2 - 7) / 9 / sqrt(2 * $n)
	- 2/405 / $n * ($u2 * (3 *$u2 + 7) - 16));
	
	if ($n <= 100) {
	my ($x0, $p1, $z);
	do {
	$x0 = $x;
	if ($x < 0) {
	$p1 = 1;
	} elsif ($n>100) {
	$p1 = _subuprob((($x / $n)**(1/3) - (1 - 2/9/$n))
	/ sqrt(2/9/$n));
	} elsif ($x>400) {
	$p1 = 0;
	} else {
	my ($i0, $a);
	if (($n % 2) != 0) {
	$p1 = 2 * _subuprob(sqrt($x));
	$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
	$i0 = 1;
	} else {
	$p1 = $a = exp(-$x/2);
	$i0 = 2;
	}
	

	for (my $i = $i0; $i <= $n-2; $i += 2) {
	$a *= $x / $i;
	$p1 += $a;
	}
	}
	$z = exp((($n-1) * log($x/$n) - log(4*PI*$x) + $n - $x - 1/$n/6) / 2);
	$x += ($p1 - $p) / $z;
	$x = sprintf("%.5f", $x);
	} while (($n < 31) && (abs($x0 - $x) > 1e-4));
	}
	}
	return $x;
}

sub log10 {
		my $n = shift;
		return log($n) / log(10);
}
 
sub max {
	my $max = shift;
	my $next;
	while (@_) {
	$next = shift;
	$max = $next if ($next > $max);
	}	
	return $max;
}

sub min {
	my $min = shift;
	my $next;
	while (@_) {
		$next = shift;
		$min = $next if ($next < $min);
	}	
	return $min;
}

sub precision {
		my ($x) = @_;
		return abs int(log10(abs $x) - SIGNIFICANT);
}

sub precision_string {
	my ($x) = @_;
		if ($x) {
		return sprintf "%." . precision($x) . "f", $x;
	} else {
		return "0";
	}
}

1;     

1;

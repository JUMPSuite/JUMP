#!/bin/env perl

use strict;
use Clone 'clone';
use Statistics::R;

our $VERSION = "1.13.001";

if (scalar(@ARGV)!=1)
{
	die "Usage: jump_i jump_i.params\n";
}

# initialization
my (%parahash);
parse_params($ARGV[0],\%parahash);

=head
# custome table: special initialization
if ( defined($parahash{customized_input}) and $parahash{customized_input} ne '0' )
{   
	print "Running customized input mode\n";
	if (defined($parahash{apply_cutoff}) and $parahash{apply_cutoff} eq '0' )
	{
		print "All filters are turned off\n";
		turnoff_filter(\%parahash);
	}
}
else 
{
	print "Running regular input mode\n";  
}
=cut

# mkdir output folder
my $outputDir="inf_$parahash{output_folder}";
if (-e $outputDir) {} else { system("mkdir $outputDir"); }
system("cp $ARGV[0] $outputDir");
my $logfile="$outputDir/jump_i.log";
open(LOG,">$logfile");

# read Data matrix
my (%data,%tmthash,@tmtarray);
print "Reading input matrix (may take a while)\n";
print LOG "Reading input matrix (may take a while)\n";
parse_input(\%parahash,\%data,\%tmthash,\@tmtarray);

# processing begins:
# matrix generation mode + heatmap
if (!defined($parahash{pretest_function}) or $parahash{pretest_function}==0)
{
	# matrix generation
	my %orig_data=%{clone(\%data)};	# back up data

	# Cut proteins/peptides by intensity: pair_cutoff_intensity
	cutIntensity(\%data,$parahash{pair_cutoff_intensity});

	# mean and SD matrix for log2ratio
	my (@meanarray,@sdarray);
	initializeArray(\@meanarray,scalar(keys %tmthash));
	initializeArray(\@sdarray,scalar(keys %tmthash));
	for ( my $i=1; $i<scalar(keys %tmthash); $i++ )
	{
		for ( my $j=$i+1; $j<=scalar(keys %tmthash); $j++ )
		{
			my ($mean,$sd)=pair_log2r_mean_SD(\%data,$i,$j,$parahash{pair_cutoff_percentage}/2);
			$meanarray[$i][$j]=sprintf "%.4f", $mean;
			$sdarray[$i][$j]=sprintf "%.4f", $sd;
		}
	}

	# printing
	print "Printing log2ratio mean matrix\n";
	print LOG "Printing log2ratio mean matrix\n";
	printMatrix(\@meanarray,"$outputDir/log2ratio_mean_matrix.txt");
	print "Printing log2ratio SD matrix\n";
	print LOG "Printing log2ratio SD matrix\n";
	printMatrix(\@sdarray,"$outputDir/log2ratio_SD_matrix.txt");
	# end of matrix generation
	# -------------------------------------------------------------------------

	# log2 ratio distribution printing
	print "Printing log2ratio distributions\n";
	print LOG "Printing log2ratio distributions\n";

	# print matrix: protein / peptides * plex
	printData(\%data, \@tmtarray, "$outputDir\/\.wholeData");
	
	# analysis in R:
	printLog2RatioDistributions($outputDir,$parahash{pair_cutoff_percentage}/2);
	
	# end of log2 ratio distribution printing
	# -------------------------------------------------------------------------

	# heatmap
	print "Printing heatmaps (may take a while)\n";
	print LOG "Printing heatmaps (may take a while)\n";
	my %data=%{clone(\%orig_data)}; # restore data

	# Cut proteins/peptides by intensity: cluster_cutoff_intensity
	cutIntensity(\%data,$parahash{cluster_cutoff_intensity});

	# Calculate RSD for each row
	my %featurehash;
	foreach my $p (keys %data) { $featurehash{$p}{rsd}=rel_std_dev(@{$data{$p}}); }

	# Calculate top % of most variable proteins / peptides
	my %rankhash; my $count=0;
	%rankhash = map {$_ => $count++} sort{$featurehash{$b}{rsd}<=>$featurehash{$a}{rsd}} keys %featurehash;
	foreach my $p (keys %featurehash) { $featurehash{$p}{rankpct}=$rankhash{$p}*100/scalar(keys %featurehash); }
	
	# cluster with top % of most variable proteins/peptides by RSD 
	foreach my $pct (sort{$b<=>$a} keys %{$parahash{cluster_top_pct}})
	{
		selectTopVariable(\%data,\%featurehash,$pct);
		#print scalar(keys %data),"\n";
		if (scalar(keys %data)>1)	# check if hash is empty (in case dataset is too small)
		{
			printHeatmap(\%data, $pct, \@tmtarray, $outputDir, $parahash{'bypass_row_clustering'});
		}
	}
	# rm 'Rplots.pdf'
	if (-e 'Rplots.pdf') { system("rm Rplots.pdf"); }
}
# test percentage: with fixed intensity
elsif ( defined($parahash{pretest_function}) and $parahash{pretest_function}==1 )
{
	# Cut proteins/peptides by intensity: pair_cutoff_intensity
	if ( !defined($parahash{pair_cutoff_intensity}) ) { die "Please specify pair_cutoff_intensity when pretest_function is set 1\n"; }
	print "Pretesting function: tesing percentage with fixed intensity ($parahash{pair_cutoff_intensity})\n";
	print LOG "Pretesting function: tesing percentage with fixed intensity ($parahash{pair_cutoff_intensity})\n";
	cutIntensity(\%data,$parahash{pair_cutoff_intensity});

	# pair analysis with different percentage cut
	my ($smp1,$smp2);
	# A) pair_min_var
	print "  For pair with min variance: $parahash{pair_min_var}\n";
	print LOG "  For pair with min variance: $parahash{pair_min_var}\n";
	($smp1,$smp2)=split(/\:/,$parahash{pair_min_var});
	pairAnalysis_testPct($smp1,$smp2);
	# B) pair_max_var
	print "  For pair with max variance: $parahash{pair_max_var}\n";
	print LOG "  For pair with max variance: $parahash{pair_max_var}\n";
	($smp1,$smp2)=split(/\:/,$parahash{pair_max_var});
	pairAnalysis_testPct($smp1,$smp2);
}
# test intensity: with fixed percentage
elsif ( defined($parahash{pretest_function}) and $parahash{pretest_function}==2 )
{
	if ( !defined($parahash{pair_cutoff_percentage}) ) { die "Please specify pair_cutoff_percentage when pretest_function is set 2\n"; }
	print "Pretesting function: tesing intensity with fixed percentage ($parahash{pair_cutoff_percentage}\%)\n";
	print LOG "Pretesting function: tesing intensity with fixed percentage ($parahash{pair_cutoff_percentage}\%)\n";

	# pair analysis with different intensity
	my ($smp1,$smp2);
	# A) pair_min_var
	print "  For pair with min variance: $parahash{pair_min_var}\n";
	print LOG "  For pair with min variance: $parahash{pair_min_var}\n";
	($smp1,$smp2)=split(/\:/,$parahash{pair_min_var});
	pairAnalysis_testInt($smp1,$smp2);
	# B) pair_max_var
	print "  For pair with max variance: $parahash{pair_max_var}\n";
	print LOG "  For pair with max variance: $parahash{pair_max_var}\n";
	($smp1,$smp2)=split(/\:/,$parahash{pair_max_var});
	pairAnalysis_testInt($smp1,$smp2);
}
else { die "Wrong setting of pretest_function!!!\n"; }
close LOG;

#---------------------------------------------------------------
sub turnoff_filter
{
	my ($parahash)=@_;

	$$parahash{remove_contaminants}=0;
	$$parahash{pair_cutoff_percentage}=0;
	$$parahash{pair_cutoff_intensity}=0;
	$$parahash{cluster_cutoff_intensity}=0;
	#$$parahash{}=;
}

sub printLog2RatioDistributions
{
	my ($outputDir,$pctCut)=@_;

	my $R = Statistics::R -> new();
	$R -> set('pctCut', $pctCut);
	$R -> set('outputDir', $outputDir);
	$R -> run(q`# load data
	dt=read.table(paste(outputDir,"/.wholeData",sep=''),head=T,sep="\t")
	m=as.matrix(dt[,2:ncol(dt)])
	rownames(m)=dt[,1]

	# cutDataEachSide
	cutDataEachSide=function(a,pctCut)
	{
	        a=sort(a)
        	cutN=round( length(a) * pctCut /100)
	        cutTmp2=(-1)*cutN
        	a=a[-1:cutTmp2]
	        cutTmp1=(length(a)-cutN+1)*(-1)
        	cutTmp2=(-1)*length(a)
	        a[cutTmp1:cutTmp2]
	}

	# print
	limit=5
	pdf(file=paste(outputDir,'/log2ratio_distributions.pdf',sep=''),width=10,height=10)
	#par(mfrow=c(10,10),mar=c(2,1,1,1))
	par(mfrow=c(ncol(m),ncol(m)),mar=c(2,1,1,1))
for (i in 1:ncol(m))
{
        for (j in 1:ncol(m))
        {
                if (i==j) # empty graph
                {
                        plot(0,0,axes = FALSE, col='white')
                        #text(0,0,colnames(m)[i],cex = 1.8)
                        text(0,0,colnames(m)[i],cex = 1.2)
                        if (i==1) # legend for cutoff
                        {
                                text(0,0.5,paste(pctCut*2,'% cutoff ->',sep=''),cex = 1, col='blue')
                        }
                        if (i==ncol(m)) # legend for cutoff
                        {
                                text(0,-0.5,'<- no cutoff',cex = 1, col='blue')
                        }
                }
                if (i<j)  # use user specified pctCut
                {

                        # Calculate log2ratio dstr
                        a=log2( ( m[,j] + 1 ) / ( m[,i] + 1 ) )
                        # Cut percentage
                        a=cutDataEachSide(a,pctCut)
                        # print density distribution
                        plot(density(a),main='',xlab='',ylab='',xlim=c(0-limit,limit))
                                                abline(v=0,col='grey',lty=2)
                }
                if (i>j)  # print all data (no cut)
                {
                        # Calculate log2ratio dstr
                        a=log2( ( m[,i] + 1 ) / ( m[,j] + 1 ) )
                        # print density distribution
                        plot(density(a),main='',xlab='',ylab='',xlim=c(0-limit,limit))
                        abline(v=0,col='grey',lty=2)

                }
        }
}
dev.off()`);
	$R -> stop();
}

sub printData
{
	my ($data,$header,$output)=@_;

	# print data to tmp file
	#open(OUT,">$outputDir\/\.wholeData");
	open(OUT,">$output");
	for (my $i=1;$i<scalar(@{$header});$i++) { print OUT "\t$$header[$i]"; }
	print OUT "\n";
	foreach my $p (keys %{$data})
	{
		print OUT "\"$p\"";
		for (my $i=1; $i<=$$data{$p}[0]; $i++) { print OUT "\t$$data{$p}[$i]"; }
		print OUT "\n";
	}
	close OUT;
}

sub printHeatmap
{
	my ($data,$pctCut,$header,$outputDir,$bypass_row_clustering)=@_;

	# print data to tmp file
	open(OUT,">$outputDir\/\.tmpdata");
	#print OUT "key";
	for (my $i=1;$i<scalar(@{$header});$i++) { print OUT "\t$$header[$i]"; }
	#print scalar(@{$header}),",$$header[10],\n";
	print OUT "\n";
	foreach my $p (keys %{$data}) 
	{ 
		print OUT "\"$p\""; 
		#print OUT "keys";  # to avoid weired simples in mod peptides
		for (my $i=1; $i<=$$data{$p}[0]; $i++) { print OUT "\t$$data{$p}[$i]"; }
		print OUT "\n";
	}
	close OUT;

	# use R to draw heatmap
	my $R = Statistics::R -> new();
	$R -> set('pctCut', $pctCut);
	$R -> set('outputDir', $outputDir);
	$R -> set('bypass_row_clustering', $bypass_row_clustering);
	$R -> run(q`data=read.table(paste(outputDir,"/.tmpdata",sep=''),head=T,sep="\t")
	m=as.matrix(data[,2:ncol(data)])
	rownames(m)=data[,1]
	library(gplots)
	par(mar=c(12,4,6,12))
	pdf(file=paste(outputDir,'/heatmap_',pctCut,'pct.pdf',sep=''),width=10,height=10)
	if (bypass_row_clustering == 1) {
		heatmap.2(log2(m),Rowv = FALSE, scale='r',main=paste('Top ',pctCut,'% most variable ',"\n",'proteins or peptides (n=',nrow(m),')',sep=''),
		density.info="n",trace="n",margins=c(12,12),labRow = NULL,
		col=rgb(
			c(seq(0,1,1/128),rep(1,129)),
			c(seq(0,1,1/128),seq(1,0,-1/128)),
			c(rep(1,129),seq(1,0,-1/128))
		)
		)
	} else {
		heatmap.2(log2(m),scale='r',main=paste('Top ',pctCut,'% most variable ',"\n",'proteins or peptides (n=',nrow(m),')',sep=''),
		density.info="n",trace="n",margins=c(12,12),labRow = NULL,
		col=rgb(
			c(seq(0,1,1/128),rep(1,129)),
			c(seq(0,1,1/128),seq(1,0,-1/128)),
			c(rep(1,129),seq(1,0,-1/128))
		)
		)
	}
	dev.off()`);
	$R -> stop();
	
}

sub selectTopVariable
{
	my ($data,$featurehash,$pct)=@_;
	foreach my $p (keys %{$data})
	{
		if ( $$featurehash{$p}{rankpct}>$pct ) { delete $$data{$p}; }
	}
}

sub printMatrix
{
	my ($array,$output)=@_;

	open(OUT,">$output");

	for ( my $i=1; $i<=scalar(@tmtarray); $i++ ) { print OUT "\t$tmtarray[$i]"; }
	print OUT "\n";
	for ( my $i=1; $i<=scalar(@tmtarray); $i++ )
	{
		print OUT "$tmtarray[$i]";
		for ( my $j=1; $j<=scalar(@{$array}); $j++ )
		{
			print OUT "\t$$array[$i][$j]";
		}
		print OUT "\n";
	}

	close OUT;
}

sub initializeArray
{
	my ($array,$n)=@_;

	for ( my $i=1; $i<=$n; $i++ )
	{
		for ( my $j=1; $j<=$n; $j++ )
		{
			$$array[$i][$j]='na';
		}
	}
}

sub pairAnalysis_testInt
{
	my ($smp1,$smp2)=@_;

	if (!defined($tmthash{$smp1}))  { die "Not defined sample: $smp1\n"; }
	if (!defined($tmthash{$smp2}))  { die "Not defined sample: $smp2\n"; }
	print "\tcutoff intensity\tn\tmean\tSD\n";
	print LOG "\tcutoff intensity\tn\tmean\tSD\n";
	for ( my $intCut=1000; $intCut<=9000; $intCut+=1000 )
	#for ( my $intCut=1000; $intCut<=90000; $intCut+=10000 )
	{
		my %tmpD=%{clone(\%data)};
		cutIntensity(\%tmpD,$intCut);
		my ($mean,$sd,$n)=pair_log2r_mean_SD(\%tmpD,$tmthash{$smp1},$tmthash{$smp2},$parahash{pair_cutoff_percentage}/2);
		printf "\t$intCut\t\t\t%d\t%.4f\t%.4f\n",$n,$mean,$sd;
		printf LOG "\t$intCut\t\t\t%d\t%.4f\t%.4f\n",$n,$mean,$sd;
	}
}

sub pairAnalysis_testPct
{
	my ($smp1,$smp2)=@_;

	if (!defined($tmthash{$smp1})) 	{ die "Not defined sample: $smp1\n"; }
	if (!defined($tmthash{$smp2})) 	{ die "Not defined sample: $smp2\n"; }
	print "\tcutoff %\tn\tmean\tSD\n";
	print LOG "\tcutoff %\tn\tmean\tSD\n";
	for ( my $pct=10; $pct<=90; $pct+=10 )
	{
		my ($mean,$sd,$n)=pair_log2r_mean_SD(\%data,$tmthash{$smp1},$tmthash{$smp2},$pct/2);
		printf "\t$pct%%\t\t%d\t%.4f\t%.4f\n",$n,$mean,$sd;
		printf LOG "\t$pct%%\t\t%d\t%.4f\t%.4f\n",$n,$mean,$sd;
		#printf "$pct%%\t\t%.4f\t%.4f\n",$mean,$sd;
	}
	
}

sub cutIntensity
{
	my ($data,$cutoff)=@_;

	foreach my $p (keys %{$data})
	{
		if ( median( @{$$data{$p}} ) < $cutoff )
		{ delete $$data{$p}; }
	}
	
}

sub pair_log2r_mean_SD
{
	my ($data,$s1,$s2,$pctCut)=@_;
	$pctCut |= 0;

	my @log2r;

	# Calculate log2ratio dstr for a specific pair
	foreach my $p (keys %{$data})
	{
		push @log2r, log( $$data{$p}[$s1] / ($$data{$p}[$s2] + 1) )/log(2);
	}

	# Cut percentage: pair_cutoff_percentage
	cutDataEachEnd(\@log2r,$pctCut);

	# return mean and SD
	return ( mean(@log2r), std_dev(@log2r), scalar(@log2r));
}

sub parse_params
{
        my ($par,$parahash)=@_;
        open(IN,$par);
        while (<IN>)
        {
                s/^\s+//;
                next if (/^#/);
                chomp;
                next unless (/ = /);

		s/\s*([;\#].*)$//; # delete comments
                my ($key,$value) = split(' = ',$_);
                if ($key =~ /^sig1[23]/)
                {
                        $$parahash{'samplelabels'}{$key}=$value;
                }
                else
                {
                        $$parahash{$key}=$value;
                }
        }
        close IN;

	# define top % variable cutoffs in %{$$parahash{cluster_top_pct}}
	#%{$$parahash{cluster_top_pct}}=('1','','5','','10','','20','');
	# custom mode
	if ( defined($$parahash{customized_input}) and $$parahash{customized_input} ne '0' )
	{  
		%{$$parahash{cluster_top_pct}}=('100','');
	}
	else # regular mdoe
	{
		%{$$parahash{cluster_top_pct}}=('0.1','','0.5','','1','','5','','10','','20','');
		$$parahash{cluster_top_pct}{$$parahash{cluster_cutoff_percentage}}='';
	}
}

sub parse_input
{
	my ($parahash,$data,$tmthash,$tmtarray)=@_;

	my $input;
	if ( defined($$parahash{customized_input}) and
	        $$parahash{customized_input} ne '0' ) # custom mode
		{ $input=$$parahash{customized_input}; }
	else { $input=$$parahash{input_table};}       # regular -q format file

	#print "Reading input file: $input\n";
	open(IN,$input) || die "Cannot open input file: $input\n";
	my $headerline=1;
	my $mode=''; my $keycol=-1; my $sigStartCol=-1; my $plex=0;
	my $con=0;
	my @dataindex2tabcol; # a temoparay hash that records: data index ~ table col
	while (<IN>)
	{
		chomp;
		chop if (/\r$/);
		#chop; # for windows ending
		if (/sig126/) # header line
		{
			$headerline=0;

			# protein or peptide mode? => defined keycol
=head
			if ( defined($$parahash{customized_input}) and 
				$$parahash{customized_input} ne '0' ) 
			{ $mode='cu'; $keycol=0; } # custom mode; check $mode
=cut
			if (/^Protein/) { $mode='pro'; $keycol=1; }
			elsif (/^Peptides/) { $mode='pep'; $keycol=0; }
			else { $mode='cu'; $keycol=0; } # custom mode; check $mode

			# define intensity col
			my @t=split(/\t/,$_);
			for (my $i=0; $i<=$#t; $i++)
			{
				#if ( $t[$i] eq 'sig126' ) { $sigStartCol=$i; last; }
				if (defined($$parahash{'samplelabels'}{$t[$i]}))
				{
					$plex++;
					$$tmthash{$$parahash{'samplelabels'}{$t[$i]}}=$plex;
					$$tmtarray[$plex]=$$parahash{'samplelabels'}{$t[$i]};

					$dataindex2tabcol[$plex]=$i;
				}
			}
=head
			# define $plex
			$plex=$#t-$sigStartCol+1;

			# define %tmthash
			my $k=1;
			for (my $i=$sigStartCol; $i<=$#t; $i++)
			{
				$$tmthash{$t[$i]}=$k; 
				$$tmtarray[$k]=$t[$i];
				$k++; #print "$t[$i]\n";
			}
=cut
		}
		elsif ( $headerline==1 ) { next; } # front lines before header
		else # real data begins
		{
			# remove_contaminants?
			if (defined($$parahash{remove_contaminants}) 
			and $$parahash{remove_contaminants} 
			and (/Contaminant/ or /contaminant/ or /co\|CON_/) )
			#{ $con++; print "$_\n";next; }
			{ $con++; next; }

			# read data
			my @t=split(/\t/,$_);
			$$data{$t[$keycol]}[0]=$plex;
			for ( my $i=1; $i<=$plex; $i++ )
			{
				#$$data{$t[$keycol]}[$i]=$t[$sigStartCol-1+$i];
				$$data{$t[$keycol]}[$i]=$t[$dataindex2tabcol[$i]];
			}
		}
	}
	close IN;

	print "Remove $con rows by contaminant filtering\n";
	print LOG "Remove $con rows by contaminant filtering\n";

	return $mode;
}

# basic statistics functions cp from http://www.perlmonks.org/?node_id=1089988
sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}
sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}
sub std_dev {
    my (@data) = @_;
    my ( $sq_dev_sum, $avg ) = ( 0, 0 );
    
    $avg = mean(@data);
    foreach my $elem (@data) {
        $sq_dev_sum += ( $avg - $elem )**2;
    }
    return ( sqrt( $sq_dev_sum / ( @data - 1 ) ) );
}
# end of cp

sub rel_std_dev  # as percentage (%)
{
	my (@data) = @_;
	return ( std_dev(@data) * 100 / mean(@data) );
}

sub cutDataEachEnd
# sort data and cut $pct (%) data at each end
{
	my ($data,$pct) = @_;

	# sort data
	@{$data} = sort {$a<=>$b} @{$data};

	# calculate how many should be cut at each end
	my $cutN=int( scalar(@{$data}) * $pct/100);

	# cut $cutN data points at each end
	splice(@{$data},0,$cutN);  # cut left/lower end
	splice(@{$data},scalar(@{$data})-$cutN); # cut right/higher end
}



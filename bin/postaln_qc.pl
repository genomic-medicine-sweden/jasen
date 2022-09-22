#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use IPC::Cmd qw[can_run run];
use JSON;# qw( encode_json );

my $VERSION = "1.0";

# DEPENDENCIES:
#  * sambamba 0.6.4+

unless( $ARGV[0] and $ARGV[1] and $ARGV[2] ) {
    print STDERR "USAGE: postaln_qc.pl BAM TARGETED_BED SAMPLE_ID [THREADS] [BAITS_BED] [REFERENCE_FASTA]\n";
    exit;
}

#die "ERROR: sambamba not in path!" unless can_run('sambamba');

my %results;
my $BAM     = $ARGV[0] or usage(1);
my $BED     = $ARGV[1] or usage(1);
my $SID     = $ARGV[2] or usage(1);
my $THREADS = ( $ARGV[3] or 0 );
my $BAITS   = ( $ARGV[4] or 0 );
my $REF_FA  = ( $ARGV[5] or 0 );
my $PAIRED  = is_PE( $BAM );


if( $BAITS and $REF_FA ) {
    print STDERR "Calculating HS-metrics...\n";
    my $DICT = $REF_FA;
    if( -s $DICT.".dict" ) {
	$DICT = $DICT.".dict";
    }
    else{
	$DICT =~ s/\.(fa|fasta)$/\.dict/;
	die "Could not find dict file for reference fasta" unless ( -s $DICT );
    }
    system_p( "picard BedToIntervalList -I $BED -O $BED.interval_list -SD $DICT" ) unless -s "$BED.interval_list";
    system_p( "picard BedToIntervalList -I $BAITS -O $BAITS.interval_list -SD $DICT" ) unless -s "$BAITS.interval_list";
    system_p( "picard CollectHsMetrics -I $BAM -O $BAM.hsmetrics -R $REF_FA -BAIT_INTERVALS $BAITS.interval_list -TARGET_INTERVALS $BED.interval_list" );

    open( HS, "$BAM.hsmetrics" );
    while( <HS> ) {
	if( /^\#\# METRICS CLASS/ ) {
	    <HS>;
	    my $vals = <HS>;
	    my @a = split /\t/, $vals;
            $results{'pct_on_target'} = $a[18];
	    $results{'fold_enrichment'} = $a[26];
	    $results{'median_coverage'} = $a[23];
	    $results{'fold_80'} = $a[33];
	}
    }
}


# Get total number of reads and number of mapped reads
print STDERR "Collecting basic stats...\n";
my $exec_str = "sambamba flagstat ".($THREADS ? "-t $THREADS": ""). " $BAM";
my @flagstat = `$exec_str`;
my( $num_reads ) = ( $flagstat[0] =~ /^(\d+)/ );
my( $dup_reads ) = ( $flagstat[3] =~ /^(\d+)/ );
my( $mapped_reads ) = ( $flagstat[4] =~ /^(\d+)/ );

if( $PAIRED ) {
    print STDERR "Collect insert sizes...\n";
    system_p( "picard CollectInsertSizeMetrics -I $BAM -O $BAM.inssize -H $BAM.ins.pdf -STOP_AFTER 1000000");
    open( INS, "$BAM.inssize" );
    while( <INS> ) {
	if( /^\#\# METRICS CLASS/ ) {
	    <INS>;
	    my $vals = <INS>;
	    my @a = split /\t/, $vals;
	    $results{'ins_size'} = $a[0];
	    $results{'ins_size_dev'} = $a[1];
	}
    }
    close INS;
    unlink( "$BAM.inssize" );
    unlink( "$BAM.ins.pdf" );
}





my $OUT_PREFIX = $BAM."_postalnQC";
my @thresholds = qw( 1 10 30 100 250 500 1000);

print STDERR "Collecting depth stats...\n";
#system_p( "sambamba depth base --fix-mate-overlaps -c 0 ".($THREADS ? "-t $THREADS": "")." -L $BED $BAM > $OUT_PREFIX.basecov.bed" );
system_p( "sambamba depth base -c 0 ".($THREADS ? "-t $THREADS": "")." -L $BED $BAM > $OUT_PREFIX.basecov.bed" );
#system_p( "samtools depth -a -b $BED ".($THREADS ? "-@ $THREADS": "")." $BAM > $OUT_PREFIX.basecov.bed" );
my( $pct_above, $mean_cov, $iqr_median ) = parse_basecov_bed( $OUT_PREFIX.".basecov.bed", \@thresholds );

unlink( $OUT_PREFIX.".basecov.bed" );


$results{pct_above_x} = $pct_above;
$results{tot_reads} = $num_reads;
$results{mapped_reads} = $mapped_reads;
$results{dup_reads} = $dup_reads;
$results{dup_pct} = $dup_reads / $mapped_reads;
$results{sample_id} = $SID;
$results{mean_cov} = $mean_cov;
$results{iqr_median} = $iqr_median;

#print encode_json(\%results);
my $json = JSON->new->allow_nonref;
print $json->pretty->encode( \%results );



sub parse_basecov_bed {
    my( $fn, $thresholds ) = @_;
    open( my $cov_fh, $fn );

    chomp( my $head_str = <$cov_fh> );
    
    $head_str =~ s/^#\s+//;
    my @head = split /\t/, $head_str;
    my $cov_field;
    for my $i ( 0..$#head ) {
	$cov_field = $i if $head[$i] eq "COV";
    }
 

    
    my $tot_bases = 0;
    my %above_cnt;
    $above_cnt{$_}=0 foreach @$thresholds;
    
    my( $tot, $cnt ) = (0,0);
    my %levels;
    while( <$cov_fh> ) {
	chomp;
	my @a = split /\t/;

	next if $a[0] =~ /^chr(Un|\d+_)/ ;

	$tot += $a[2];
	$cnt++;

	$levels{ $a[ $cov_field ] } ++;
	
	$tot_bases ++;
	foreach my $min ( @$thresholds ) {
	    $above_cnt{ $min } ++ if $a[ $cov_field ] >= $min; 
	}
    }

    my %above_pct;
    #print "Total: $tot_bases\n";
    foreach( sort {$a<=>$b} keys %above_cnt ) {
	$above_pct{ $_ } = 100 * ($above_cnt{$_} / $tot_bases);
    }

    my $mean_cov = $tot / $cnt;

    # Calculate the inter-quartile range / median (IQR/median)
    my $q1_num = $cnt / 4;
    my $q3_num = 3 * $cnt / 4;
    my $median_num = $cnt / 2;
    my $sum = 0;
    my( $q1, $q3, $median );
    my $iqr_median = "9999";
    foreach my $l ( sort {$a<=>$b} keys %levels ) {
	$sum += $levels{$l};
	
	if( $sum >= $q1_num and !$q1 ) {
	    $q1 = $l;
	}
	if( $sum >= $median_num and !$median ) {
	    $median = $l;
	}
	if( $sum >= $q3_num and !$q3 ) {
	    $q3 = $l;
	}
    }
    if( $q1 and $q3 and $median ) {
	$iqr_median = ($q3-$q1) / $median;
    }
    
    return \%above_pct, $mean_cov, $iqr_median;
}



sub parse_cov_bed {
    my( $fn, $thresholds ) = @_;
    open( my $cov_fh, $fn );

    chomp( my $head_str = <$cov_fh> );
    
    $head_str =~ s/^#\s+//;
    my @head = split /\t/, $head_str;
    my $cov_field;
    for my $i ( 0..$#head ) {
	$cov_field = $i if $head[$i] eq "meanCoverage";
    }
 

    my $tot_exons = 0;
    my %above_cnt;
    while( <$cov_fh> ) {
	chomp;
	my @a = split /\t/;
	$tot_exons ++;
	foreach my $min ( @$thresholds ) {
	    $above_cnt{ $min } ++ if $a[ $cov_field ] >= $min; 
	}
    }

    foreach( sort {$a<=>$b} keys %above_cnt ) {
	print "$_\t";
	printf "%.2f%%", 100 *($above_cnt{$_} / $tot_exons);
	print "\n";
    }

}






sub usage {
    my $kill = shift;
    print "postaln_qc.pl BAM BED SE/PE [threads] \n";
    exit if $kill;
}


sub system_p {
    my @cmd = @_;
    print STDERR "RUNNING: ".join " ", @cmd;
    print STDERR "\n";
    system( join " ", @cmd );
}


sub is_PE {
    my $bamfile = shift;

    chomp( my $line=`samtools view $bamfile | head -n 1| awk '{print \$2}'`);
    my $remainder=$line%2;

    my $is_paired = ( $remainder ? 1 : 0 );

    return $is_paired;
}


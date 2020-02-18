#!/usr/bin/perl -w
use strict;
use JSON;
use Data::Dumper;
use CMD::tsv qw(read_tsv);

my $db_fn = $ARGV[0];
my $summary_fn = $ARGV[1];
my $report_fn = $ARGV[2];
my %results;


# Get genes that are in the database into a hash
open( DB, $db_fn );
while(<DB>) {
    chomp;
    if( /^>(.*?)\./ ) {
	$results{$1} = {'present' => 0}
    }
}
close DB;

# Get summary information (present/absent) and update the hash
open( SUM, $summary_fn );
chomp( my $head = <SUM> );
chomp( my $data = <SUM> );
my @head = split /,/, $head;
my @data = split /,/, $data;

for(1..$#head) {
    my $id = $head[$_];
    $id =~ s/\.match//;
    if( $data[$_] eq "yes" ) {
	$results{$id}->{'present'} = 1;
    }
}
close SUM;


# Add identity% and gene coverage data
my @report = read_tsv($report_fn);
foreach my $row ( @report ) {
    my $gene   = $row->{cluster};
    my $contig = $row->{ctg};
    $contig =~ s/\./_/g;
    die "Gene in report that does not exist in database: $gene" unless $results{$gene};

    unless( $results{$gene}->{$contig} ) {
	$results{$gene}->{$contig}->{id} = $row->{pc_ident};
	$results{$gene}->{$contig}->{ref_len} = $row->{ref_len};
	$results{$gene}->{$contig}->{match_len} = $row->{ref_base_assembled};
    }
}


# Print to JSON
my $json = JSON->new->allow_nonref;
print $json->pretty->encode( \%results );

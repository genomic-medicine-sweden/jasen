#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $fasta = $ARGV[0];
my $vcf = $ARGV[1];

my %fa = read_fasta($fasta);
my %sites = read_vcf($vcf);

mask_sites(\%fa, \%sites);

foreach ( keys %fa ) {
    print ">$_\n".$fa{$_}."\n";
}

sub mask_sites {
    my ($fa, $sites) = @_;

    foreach ( keys %$sites ) {
	my( $contig, $pos ) = split /:/;
	substr( $fa->{$contig}, $pos-1, 1, "N" );
    }
}
sub read_fasta {
    my $fn = shift;
    open( FA, $fn );
    my %fa;
    my $id;
    while(<FA>) {
	chomp;
	next if /^\s*$/;
	if( /^>(.*?)$/ ) {
	    $id = $1;
	}
	else {
	    $fa{$id} .= $_;
	}
    }
    close FA;
    return %fa;
}

sub read_vcf {
    my $fn = shift;
    open( VCF, $fn );
    my %sites;
    while(<VCF>) {
	chomp;
	next if /^#/;
	my @a = split /\t/;
	$sites{$a[0].":".$a[1]} = 1;
    }
    close VCF;
    return %sites;
}

#!/usr/bin/perl -w
use strict;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;
use Data::Dumper;
use Getopt::Long;
use JSON;
use Digest::MD5 qw( md5_hex );
use File::Basename qw( basename );
use CMD::tsv qw(read_tsv);

# TODO: Implement --overwrite in a reasonable way...

my %opt;
GetOptions( \%opt, 'in=s', 'id=s', 'run=s', 'overwrite', 'metadata=s', 'species=s', 'mlst=s', 'kraken=s', 'aribavir=s', 'quast=s', 'dry' );


print_usage(1) unless ( $opt{in} or $opt{metadata} ) and $opt{run} and $opt{species};

die "File not found: $opt{in}\n" if $opt{in} and ! -e $opt{in};
die "File not found: $opt{metadata}\n" if $opt{metadata} and ! -e $opt{metadata};
die "File not found: $opt{mlst}\n" if $opt{mlst} and ! -e $opt{mlst};
die "File not found: $opt{kraken}\n" if $opt{kraken} and ! -e $opt{kraken};
die "File not found: $opt{aribavir}\n" if $opt{aribavir} and ! -e $opt{aribavir};

# Connet to database, and create handles for collections
my $client = MongoDB->connect();
my $SCHEME = $client->ns("cgviz.scheme");
my $SAMPLE = $client->ns("cgviz.sample");


if( $opt{in} ) {
    my ( $scheme_md5, $scheme, $data, $missing ) = read_chewbbaca( $opt{in} );

    # Read metadata if a file was given
    my $metadata;
    if( $opt{metadata} ) {
	$metadata = read_metadata( $opt{metadata} );
    }

    # Read MLST data if file was given
    my $mlst;
    if( $opt{mlst} ) {
	$mlst = read_mlst( $opt{mlst} );
    }

    # Read brakken (kraken) data if file was given
    my( $brakken, $top_brakken );
    if( $opt{kraken} ) {
	($brakken, $top_brakken) = read_brakken( $opt{kraken} );
    }

    # Read ariba virulence data if file was given
    my $aribavir;
    if( $opt{aribavir} ) {
	$aribavir = read_aribavir( $opt{aribavir} );
    }

    # Read quast assembly QC data if file was given
    my $quast;
    if( $opt{quast} ) {
	$quast = read_quast( $opt{quast} );
    }
    
    # Create new scheme if the one from the input file doesn't already exist
    unless( scheme_exists($scheme_md5, $opt{species}) ) {
	print "New CGMLST scheme! Adding to database.\n";
	$SCHEME->insert_one( {'md5'=>$scheme_md5, 'loci'=>$scheme, 'species'=>$opt{species} } ) unless $opt{dry};
    }

    # Insert allele data for single samples
    if( scalar keys %$data == 1 ) {
	my $id = ( $opt{id} or clean_id( (keys %$data)[0] ) );
	my $exists = sample_exists( $scheme_md5, $opt{species}, $id, $opt{run} );

	my $run_name = (split '/', $opt{run})[-1];

	my $run_date = (split '_', $run_name)[-1];
	
	# If sample already exists and overwrite was requested, completely delete entry from database
	if( $exists and $opt{overwrite} ) {
	    $SAMPLE->delete_one( {'scheme_md5'=>$scheme_md5, 'id'=>$id, 'species'=>$opt{species}, 'run'=>$opt{run}} ) unless $opt{dry};
	    $exists = 0;
	}

	if( !$exists ) { 
	    $SAMPLE->insert_one(
		{
		    'scheme_md5' => $scheme_md5,
		    'alleles'    => (values %$data)[0], 
		    'id'         => $id, 
		    'species'    => $opt{species}, 
		    'run'        => $opt{run}, 
		    'metadata'   => ( $metadata->{$id} or {} ), 
		    'missing'    => (values %$missing)[0], 
		    'mlst'       => $mlst,
		    'brakken'    => $brakken,
		    'top_brakken'=> $top_brakken, 
		    'quast'      => $quast,
		    'aribavir'   => $aribavir
		} ) unless $opt{dry};
	}
	else {
	    print "Sample with same sampleID, runID, species and cgmlst scheme already exists in database ($id). Skipping! Use --overwrite to override\n";
	}
    }
    # Insert allele data for multiple samples
    else {
	die "Multiple samples in chewbbaca file! Use import_cgmlst_multi.pl...";
    }
}




sub clean_id {
    my $id = shift;
    $id =~ s/\.spades//;
    return $id;
}

sub scheme_exists {
    my $md5 = shift;
    my $species = shift;
    my $results = $SCHEME->find_one( {'md5'=>$md5, 'species'=>$species } );
    return 1 if $results and $results->{md5} eq $md5;
    return 0;  
}

sub sample_exists {
    my $md5     = shift;
    my $species = shift;
    my $id      = shift;
    my $run     = shift;
    my $results = $SAMPLE->find_one( {'scheme_md5'=>$md5, 'species'=>$species,  'id'=>$id, 'run'=>$opt{run} } );
    return 1 if $results and $results->{id} eq $id;
    return 0;  
}

sub read_chewbbaca {
    my $fn = shift;

    open( CHEW, $fn ) or die "Cannot open file: $fn\n";
    chomp( my $head = <CHEW> );
    my $head_md5 = md5_hex($head.$opt{species});
    my @head = split /\t/, $head;
    shift @head;

    my %result;
    my %missing;
    while( <CHEW> ) {
	chomp;
	my @data = split /\t/;
	my $id = shift @data;
	die "Error: Chewbbaca file contains same sample multiple times ($id)" if $result{$id};
        $missing{clean_id($id)} = grep /^-$/, @data;
	$result{$id} = \@data;
    }
    return $head_md5, \@head, \%result, \%missing;
}

sub read_metadata {
    my $fn = shift;

    open( META, $fn ) or die "Cannot open file: $fn\n";
    chomp( my $head = <META> );
    my @head = split /\t/, $head;

    my %metadata;
    while( <META> ) {
	chomp;
	my @data = split /\t/;
	my $id = $data[0];
	my %row;
	for my $i ( 1..$#data ) {
	    $row{$head[$i]} = $data[$i];
	}
	$metadata{$id} = \%row;
    }
    return \%metadata;
}


sub read_mlst {
    my $fn = shift;
    my $json = read_json($fn);
    my $data;
    foreach my $key (keys %{$json}) {
	my $s = $json->{$key};
	$data = {'sequence_type' => $s->{sequence_type},
		 'alleles'       => $s->{alleles} } 
    }
    return $data;
}


sub read_brakken {
    my $fn = shift;
    my @brakken = read_tsv($fn);

    my $top = (sort {$b->{fraction_total_reads} <=> $a->{fraction_total_reads}} @brakken)[0];
    
    return \@brakken, {'top_species'=>$top->{name}, 'top_species_pct'=>$top->{fraction_total_reads}};
}

sub read_quast {
    my $fn = shift;
    my @quast = read_tsv($fn);

    # Get rid of the "." in on of the keys, which is not allowed by mongodb.
    $quast[0]->{'# unaligned mis contigs'} = $quast[0]->{'# unaligned mis. contigs'};
    delete $quast[0]->{'# unaligned mis. contigs'};
    
    return $quast[0];    
}


sub read_aribavir {
    my $fn = shift;
    my $json = read_json($fn);
    return $json;
}

sub print_usage {

    print "\nUSAGE: import_cgmlst.pl\n\n".
	  "  --in TSV         CGMLST profiles from chewbbaca (required)\n".
	  "  --metadata TSV   File with metadata\n".
	  "  --mlst JSON      File with 7-gene MLST data\n".	  
	  "  --id ID          Sample ID (ID from chewbbaca file used if not given)\n".
	  "  --run ID         Run ID (required)\n".
	  "  --species ID     Species name (required)\n".
	  "  --dry            Don't write anything to the database\n".
	  "  --overwrite      Overwrite data in database (default: do not overwrite) ONLY WORKS FOR SINGLE SAMPLE FILES!\n\n";
    
    exit if $_[0];
}



sub read_json {
    my $fn = shift;

    print STDERR "Reading json $fn\n";

    open( JSON, $fn );
    my @json = <JSON>;
    my $decoded = decode_json( join("", @json ) );
    close JSON;

    return $decoded;
}

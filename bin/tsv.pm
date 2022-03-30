package tsv;
use strict;
use Exporter qw(import);

our @EXPORT_OK = qw( read_tsv ); 



# Parse VCF file and return
sub read_tsv {
    my $fn = $_[0];
    die "File does not exist: $fn" unless -e $fn;

    my $outfmt = "array";
    my $key;
    if( $_[1] and $_[1] ne "_ARRAY" ) {
	$key = $_[1];
	$outfmt = "hash";
    }

    my $lenient = 0;
    $lenient = 1 if $_[2] and $_[2] eq "lenient";
    
    open( TSV, $fn );

    chomp( my $header = <TSV> );
    $header =~ s/\r//g;
    my @head = split /\t/, $header;

    # Make sure that specified key field exists in header, and save position
    my $key_column_num;
    if( $outfmt eq "hash" ) {
	for( 0..$#head ) {
	    $key_column_num = $_ if $head[$_] eq $key;
	}
	die "The header field '$key' does not exist (".join(", ", @head).")" unless defined $key_column_num;
    }

    
    my %out_hash;
    my @out_arr;
    while( <TSV> ) {
	chomp;

	s/\r//g;
	
	my @data = split /\t/;

	# Check if equal number of columns as header
	die "Data row does not have the same number of columns as header..." unless $#data == $#head or $lenient;

	# Save row data to row hash
	my %row;
	for my $i ( 0..$#data ) {
	    $row{ $head[$i] } = $data[$i];
	}

	if( $outfmt eq "hash" ) {
	    die "Key value is not unique: ".$data[$key_column_num] if $out_hash{ $data[$key_column_num] };
	    $out_hash{ $data[$key_column_num] } = \%row;
	    
	}
	else {
	    push( @out_arr, \%row );
	}
    }


    return %out_hash if $outfmt eq "hash";
    return @out_arr if $outfmt eq "array";
    
}

1;

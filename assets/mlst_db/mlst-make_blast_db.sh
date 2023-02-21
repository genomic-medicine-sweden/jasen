#!/bin/bash

MLSTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PUBMLSTDIR="$MLSTDIR/pubmlst"
BLASTDIR="$MLSTDIR/blast"
BLASTFILE="$BLASTDIR/mlst.fa"

mkdir -p "$BLASTDIR"
rm -f "$BLASTFILE"

for N in $(find $PUBMLSTDIR -mindepth 1 -maxdepth 1 -type d); do
  SCHEME=$(basename $N)
  echo "Adding: $SCHEME"
  cat "$PUBMLSTDIR"/$SCHEME/*.tfa \
  	| grep -v 'not a locus'  \
  	| sed -e "s/^>/>$SCHEME./" \
  	>> "$BLASTFILE"
done

makeblastdb -hash_index -in "$BLASTFILE" -dbtype nucl -title "PubMLST" -parse_seqids

echo "Created BLAST database for $BLASTFILE"


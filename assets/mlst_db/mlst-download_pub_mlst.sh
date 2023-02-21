#!/bin/bash

set -e

MLSTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR="$MLSTDIR/pubmlst"
mkdir -p "$OUTDIR"
wget --no-clobber -P "$OUTDIR" http://pubmlst.org/data/dbases.xml

for URL in $(grep '<url>' $OUTDIR/dbases.xml); do
#  echo $URL
  URL=${URL//<url>}
  URL=${URL//<\/url>}
#  echo ${URL: -4}
  if [ ${URL:(-4)} = "_csv" ]; then
    #PROFILE=$(basename $URL .txt)
    PROFILE=$(echo $URL | awk -F'_' '{print $2}')
    if [ $(echo $URL | awk -F'/' '{print $3}')  = "rest.pubmlst.org" ]; then
        NUM=$(echo $URL | awk -F'/' '{if($7!=1) print "_"$7}')
    else
        NUM=$(echo $URL | awk -F'/' '{if($8!=1) print "_"$8}')
    fi
    echo "# $PROFILE "
    PROFILEDIR="$OUTDIR/$PROFILE$NUM"
    eval "mkdir -p '$PROFILEDIR'"
    eval "(cd '$PROFILEDIR' && echo "$URL" && wget -q '$URL' -O '$PROFILE$NUM.txt')"
  elif [ ${URL:(-6)} = "_fasta" ]; then
    if [ $(echo $URL | awk -F'/' '{print $3}')  = "rest.pubmlst.org" ]; then
        ALLELE=$(echo $URL | awk -F'/' '{print $7}')
    else
        ALLELE=$(echo $URL | awk -F'/' '{print $8}')
    fi
    eval "(cd '$PROFILEDIR' && echo "$URL" && wget -q '$URL' -O '$ALLELE.tfa')"
  fi
done

# delete fungi schemes
echo rm -frv "$OUTDIR"/{afumigatus,blastocystis,calbicans,cglabrata,ckrusei}
echo rm -frv "$OUTDIR"/{ctropicalis,csinensis,kseptempunctata,sparasitica,tvaginalis}
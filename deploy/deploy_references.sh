#!/bin/bash

mkdir assets &> /dev/null
mkdir assets/genomes &> /dev/null
mkdir assets/card &> /dev/null
mkdir assets/cgmlst &> /dev/null
mkdir assets/blast &> /dev/null

source activate jasen

cd assets/card
wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xjf broadstreet-v3.1.4.tar.bz2
#wget https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2
cd ../..
cd assets/blast
wget https://raw.githubusercontent.com/tseemann/mlst/master/db/blast/mlst.fa --no-check-certificate
cd ../..

#SAureus
python bin/download_ncbi.py NC_007795.1 assets/genomes
mkdir -p assets/cgmlst/staphylococcus_aureus/alleles &> /dev/null
cd assets/cgmlst/staphylococcus_aureus/alleles  
wget https://www.cgmlst.org/ncs/schema/141106/alleles/ --no-check-certificate
unzip index.html
cd ../../../..

#EColi
python bin/download_ncbi.py GCF_000008865.2 assets/genomes
mkdir -p assets/cgmlst/escherichia_coli/alleles &> /dev/null
cd assets/cgmlst/escherichia_coli/alleles
wget https://www.cgmlst.org/ncs/schema/5064703/alleles/ --no-check-certificate
unzip index.html
cd ../../../..

#KPneumoniae
python bin/download_ncbi.py NC_016845.1 assets/genomes
mkdir -p assets/cgmlst/klebsiella_pneumoniae/alleles &> /dev/null
cd assets/cgmlst/klebsiella_pneumoniae/alleles
wget https://www.cgmlst.org/ncs/schema/2187931/alleles/ --no-check-certificate
unzip index.html
cd ../../../..

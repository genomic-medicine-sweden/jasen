#!/bin/bash

mkdir assets &> /dev/null
mkdir assets/genomes &> /dev/null
mkdir assets/card &> /dev/null
mkdir assets/cgmlst &> /dev/null

source activate jasen

cd assets/card
wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xjf broadstreet-v3.1.4.tar.bz2
#wget https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2
cd ../..

#SAureus
python bin/download_ncbi.py NC_007795.1 assets/genomes
mkdir assets/cgmlst/stapylococcus_aureus &> /dev/null
cd assets/cgmlst/stapylococcus_aureus  
wget https://www.cgmlst.org/ncs/schema/141106/alleles/
unzip index.html
cd ../..

#EColi
python bin/download_ncbi.py GCF_000008865.2 assets/genomes
mkdir assets/cgmlst/escherichia_coli &> /dev/null
cd assets/cgmlst/escherichia_coli
wget https://www.cgmlst.org/ncs/schema/5064703/alleles/
unzip index.html
cd ../..

#KPneumoniae
python bin/download_ncbi.py NC_016845.1 assets/genomes
mkdir assets/cgmlst/klebsiella_pneumoniae &> /dev/null
cd assets/cgmlst/klebsiella_pneumoniae
wget https://www.cgmlst.org/ncs/schema/2187931/alleles/
unzip index.html
cd ../..

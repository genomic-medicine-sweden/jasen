#!/bin/bash

mkdir assets &> /dev/null
mkdir assets/genomes &> /dev/null
mkdir assets/card &> /dev/null
mkdir assets/cgmlst &> /dev/null
mkdir assets/blast &> /dev/null

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
assdir="${scriptdir}/../assets/"

source activate jasen

#Card
cd assets/card
wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xjf broadstreet-v3.1.4.tar.bz2
ariba prepareref -f nucleotide_fasta_protein_homolog_model.fasta --all_coding yes --force tmpdir
cp tmpdir/* .

#Mlst
#wget https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2
cd ../..
cd assets/blast
wget https://raw.githubusercontent.com/tseemann/mlst/master/db/blast/mlst.fa --no-check-certificate
cd ../..

#SAureus
python bin/download_ncbi.py CP000046.1 assets/genomes
cd assets/genomes
bwa index CP000046.1
cd ../..
mkdir -p assets/cgmlst/staphylococcus_aureus/alleles &> /dev/null
cd assets/cgmlst/staphylococcus_aureus/alleles  
wget https://www.cgmlst.org/ncs/schema/141106/alleles/ --no-check-certificate
unzip index.html
cd ..
echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee"
chewBBACA.py PrepExternalSchema  -i ${assdir}/cgmlst/staphylococcus_aureus/alleles -o ${assdir}/cgmlst/staphylococcus_aureus/alleles_rereffed \
	--cpu 1 --ptf ${assdir}/prodigal_training_files/Staphylococcus_aureus.trn
cd ../../..

#EColi
python bin/download_ncbi.py NC_000913.3 assets/genomes
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

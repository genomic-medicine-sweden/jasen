#!/bin/bash

mkdir assets &> /dev/null
mkdir assets/genomes &> /dev/null
mkdir assets/card &> /dev/null
mkdir assets/cgmlst &> /dev/null
mkdir assets/blast &> /dev/null

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
assdir="${scriptdir}/../assets/"

source activate jasen

## DBS

#CARD db
cd ${assdir}/card
wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xjf broadstreet-v3.1.4.tar.bz2
ariba prepareref -f nucleotide_fasta_protein_homolog_model.fasta --all_coding yes --force tmpdir
cp tmpdir/* .

#MLST db
#wget https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2
cd ${assdir}/blast
wget https://raw.githubusercontent.com/tseemann/mlst/master/db/blast/mlst.fa --no-check-certificate

#Finder dbs
cd ${assdir}/kma && make
cd ${assdir}/virulencefinder_db
export PATH=$PATH:/${assdir}/kma
python INSTALL.py ${assdir}/kma/kma_index
cd ${assdir}/resfinder_db
python INSTALL.py ${assdir}/kma/kma_index
cd ${assdir}/pointfinder_db
python INSTALL.py ${assdir}/kma/kma_index

## Organisms

#SAureus
cd ${assdir}/..
python bin/download_ncbi.py CP000046.1 assets/genomes
cd ${assdir}/genomes
bwa index CP000046.1.fasta
mkdir -p ${assdir}/cgmlst/staphylococcus_aureus/alleles &> /dev/null
cd ${assdir}/cgmlst/staphylococcus_aureus/alleles  
wget https://www.cgmlst.org/ncs/schema/141106/alleles/ --no-check-certificate
unzip index.html
cd ${assdir}/cgmlst/staphylococcus_aureus/ 
echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee"
chewBBACA.py PrepExternalSchema  -i ${assdir}/cgmlst/staphylococcus_aureus/alleles -o ${assdir}/cgmlst/staphylococcus_aureus/alleles_rereffed \
	--cpu 1 --ptf ${assdir}/prodigal_training_files/Staphylococcus_aureus.trn

#EColi
cd ${assdir}/..
python bin/download_ncbi.py NC_000913.3 assets/genomes
mkdir -p assets/cgmlst/escherichia_coli/alleles &> /dev/null
cd assets/cgmlst/escherichia_coli/alleles
wget https://www.cgmlst.org/ncs/schema/5064703/alleles/ --no-check-certificate
unzip index.html

#KPneumoniae
cd ${assdir}/..
python bin/download_ncbi.py NC_016845.1 assets/genomes
mkdir -p assets/cgmlst/klebsiella_pneumoniae/alleles &> /dev/null
cd assets/cgmlst/klebsiella_pneumoniae/alleles
wget https://www.cgmlst.org/ncs/schema/2187931/alleles/ --no-check-certificate
unzip index.html

cd ${assdir}/..

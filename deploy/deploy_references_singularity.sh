#!/bin/bash

mkdir assets &> /dev/null
mkdir -p assets/genomes/{escherichia_coli,klebsiella_pneumoniae,staphylococcus_aureus} &> /dev/null
mkdir assets/card &> /dev/null
mkdir assets/cgmlst &> /dev/null
mkdir -p assets/amrfinder_db/allele_counts_by_year &> /dev/null
mkdir -p assets/mlst_db/{blast,pubmlst} &> /dev/null

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
assdir="${scriptdir}/../assets/"

containerdir="${scriptdir}/../container/"
mntroot=/$(readlink -f .| cut -d"/" -f2)

## DBS

#AMR
singularity exec --bind $mntroot ${containerdir}/amrfinderplus.sif amrfinder_update -d ${assdir}/amrfinder_db

#MLST db
cd ${assdir}/mlst_db
bash ./mlst-download_pub_mlst.sh &> /dev/null
bash ./mlst-make_blast_db.sh &> /dev/null

#Finder dbs
cd ${assdir}/kma && make
cd ${assdir}/virulencefinder_db
export PATH=$PATH:/${assdir}/kma
python3 INSTALL.py ${assdir}/kma/kma_index
cd ${assdir}/resfinder_db
python3 INSTALL.py ${assdir}/kma/kma_index
cd ${assdir}/pointfinder_db
python3 INSTALL.py ${assdir}/kma/kma_index

## Organisms

#Saureus
cd ${assdir}/..
python3 bin/download_ncbi.py NC_002951.2 assets/genomes/staphylococcus_aureus
cd ${assdir}/genomes/staphylococcus_aureus
singularity exec --bind $mntroot ${containerdir}/bwakit.sif bwa index NC_002951.2.fasta
mkdir -p ${assdir}/cgmlst/staphylococcus_aureus/alleles &> /dev/null
cd ${assdir}/cgmlst/staphylococcus_aureus/alleles  
wget https://www.cgmlst.org/ncs/schema/141106/alleles/ --no-check-certificate &> /dev/null
unzip index.html &> /dev/null
cd ${assdir}/cgmlst/staphylococcus_aureus/ 
echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee"
singularity exec --bind $mntroot ${containerdir}/chewbbaca.sif chewie PrepExternalSchema -i ${assdir}/cgmlst/staphylococcus_aureus/alleles -o ${assdir}/cgmlst/staphylococcus_aureus/alleles_rereffed \
	--cpu 1 --ptf ${assdir}/prodigal_training_files/Staphylococcus_aureus.trn

#Ecoli
cd ${assdir}/..
python3 bin/download_ncbi.py NC_000913.3 assets/genomes/escherichia_coli
cd ${assdir}/genomes/escherichia_coli
singularity exec --bind $mntroot ${containerdir}/bwakit.sif bwa index NC_000913.3.fasta
mkdir -p ${assdir}/cgmlst/escherichia_coli/alleles &> /dev/null
cd ${assdir}/cgmlst/escherichia_coli/alleles
wget https://www.cgmlst.org/ncs/schema/5064703/alleles/ --no-check-certificate &> /dev/null
unzip index.html &> /dev/null
cd ${assdir}/cgmlst/escherichia_coli/
echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee"
singularity exec --bind $mntroot ${containerdir}/chewbbaca.sif chewie PrepExternalSchema -i ${assdir}/cgmlst/escherichia_coli/alleles -o ${assdir}/cgmlst/escherichia_coli/alleles_rereffed \
	--cpu 1 --ptf ${assdir}/prodigal_training_files/Escherichia_coli.trn

#Kpneumoniae
cd ${assdir}/..
python3 bin/download_ncbi.py NC_016845.1 assets/genomes/klebsiella_pneumoniae
cd ${assdir}/genomes/klebsiella_pneumoniae
singularity exec --bind $mntroot ${containerdir}/bwakit.sif bwa index NC_016845.1.fasta
mkdir -p ${assdir}/cgmlst/klebsiella_pneumoniae/alleles &> /dev/null
cd ${assdir}/cgmlst/klebsiella_pneumoniae/alleles
wget https://www.cgmlst.org/ncs/schema/2187931/alleles/ --no-check-certificate &> /dev/null
unzip index.html &> /dev/null

cd ${assdir}/..

#chewbbaca check
saureus=${assdir}/cgmlst/staphylococcus_aureus/alleles_rereffed
if [ -d "$saureus" ]; then echo "$saureus exists."; else echo "ERROR: $saureus does not exist!!! Please report this to JASEN issues."; fi

#bwa check
ref=${assdir}/genomes/staphylococcus_aureus/NC_002951.2.fasta; refamb=$ref.amb; refann=$ref.ann; refbwt=$ref.bwt; refpac=$ref.pac; refsa=$ref.sa
if [[ -f $ref && -f $refamb && -f $refann && -f $refbwt && -f $refpac && -f $refsa ]]; then echo "bwa indexes exists."; else echo "ERROR: bwa indexes do not exist!!! Please report this to JASEN issues."; fi

#blastdb check
mlst=${assdir}/mlst_db/blast/mlst.fa; mlstndb=$mlst.ndb; mlstnhd=$mlst.nhd; mlstnhi=$mlst.nhi; mlstnhr=$mlst.nhr; mlstnin=$mlst.nin; mlstnog=$mlst.nog; mlstnos=$mlst.nos; mlstnot=$mlst.not; mlstnsq=$mlst.nsq; mlstntf=$mlst.ntf; mlstnto=$mlst.nto
if [[ -f $mlst && -f $mlstndb && -f $mlstnhd && -f $mlstnhi && -f $mlstnhr && -f $mlstnin && -f $mlstnog && -f $mlstnos && -f $mlstnot && -f $mlstnsq && -f $mlstntf && -f $mlstnto ]]; then echo "BLAST indexes exists!"; else echo "ERROR: BLAST indexes do not exist!!! Please report this to JASEN issues."; fi

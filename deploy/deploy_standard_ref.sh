mkdir ../assets
mkdir ../assets/genomes
mkdir ../assets/card
mkdir ../assets/cgmlst

bash ../bin/download_ncbi.py NC_007795.1 ../assets/genomes

wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xvzf broadstreet-v3.1.4.tar.bz2
mv card-data ../assets/card
#wget https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2

wget https://www.cgmlst.org/ncs/schema/141106/alleles/
unzip Staphylococcus_aureus_cgMLST_alleles.zip
mv Staphylococcus_aureus_cgMLST_alleles ../assets/cgmlst/stapylococcus_aureus

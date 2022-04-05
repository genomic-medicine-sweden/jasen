mkdir assets
mkdir assets/genomes
mkdir assets/card
mkdir assets/cgmlst

bash bin/download_ncbi.py NC_007795.1 assets/genomes

wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xvzf broadstreet-v3.1.4.tar.bz2
mv card-data assets/card
#wget https://card.mcmaster.ca/download/5/ontology-v3.1.4.tar.bz2

#SAureus
bash bin/download_ncbi.py NC_007795.1 assets/genomes
mkdir assets/cgmlst/stapylococcus_aureus
wget https://www.cgmlst.org/ncs/schema/141106/alleles/
unzip Staphylococcus_aureus_cgMLST_alleles.zip
mv Staphylococcus_aureus_cgMLST_alleles assets/cgmlst/stapylococcus_aureus/alleles

#EColi
bash bin/download_ncbi.py GCF_000008865.2 assets/genomes
mkdir assets/cgmlst/escherichia_coli
wget https://www.cgmlst.org/ncs/schema/5064703/alleles/
unzip Escherichia_coli_cgMLST_alleles.zip
mv Escherichia_coli_cgMLST_alleles assets/cgmlst/escherichia_coli/alleles

#KPneumoniae
bash bin/download_ncbi.py NC_016845.1 assets/genomes
mkdir assets/cgmlst/klebsiella_pneumoniae
wget https://www.cgmlst.org/ncs/schema/2187931/alleles/
unzip Klebsiella_pneumoniae_cgMLST_alleles.zip
mv Klebsiella_pneumoniae_cgMLST_alleles assets/cgmlst/klebsiella_pneumoniae/alleles



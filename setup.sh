set -e
shopt -s nullglob

#Suggests provided branch. Else suggests master
default_branch=${1-master}
default_name=${2-d_jasen}

echo "Welcome to the dependency installation script. Q to exit"
while true; do
    echo "Name your conda environment ['d_jasen']:"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input = "y" ]] || [[ $input = "yes" ]]; then
        cname=$default_name
        break
    else
        cname=$input
        break
    fi
done
echo "Thank you, setting up environment $cname!"

#Unload environment
conda info| grep -q $cname && source deactivate || :
#Remove environment if already present
conda remove -y -n $cname --all || :

conda create -y -n $cname python=3.6
source activate $cname
conda config --add channels bioconda conda-forge
conda install -y -c conda-forge r-base
conda install -y -c bioconda blast=2.9.0 bwa=0.7.17 fastqc mlst mummer kraken2 picard=2.20.3 pigz=2.4 quast=5.0.2 samtools=1.9=h8571acd_11 spades=3.13.1 trimmomatic=0.39
pip install biopython ariba multiqc pandas
echo "Dependency Setup Complete!"

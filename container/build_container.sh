#!/usr/bin/env bash
set -e

definitions=(chewbbaca bedtools sambamba postAlignQc resfinder virulencefinder)
declare -A containers=( 
    [bwa]=https://depot.galaxyproject.org/singularity/bwa:0.7.17--pl5.22.0_2
    [kraken]=https://depot.galaxyproject.org/singularity/kraken:1.1.1--pl5262h7d875b9_5
    [kraken2]=https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5262h7d875b9_0
    [bracken]=https://depot.galaxyproject.org/singularity/bracken:2.6.1--py39h7cff6ad_2
    [samtools]=https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0
    [ariba]=https://depot.galaxyproject.org/singularity/ariba:2.14.6--py38h6ed170a_0
    [mlst]=https://depot.galaxyproject.org/singularity/mlst:2.19.0--hdfd78af_1
    [spades]=https://depot.galaxyproject.org/singularity/spades:3.15.2--h95f258a_1
    [quast]=https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl5262h190e900_4
    [freebayes]=https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py39hba5d119_3
    [perl_json]=https://depot.galaxyproject.org/singularity/perl-json%3A4.02--pl526_0
)

for tool in "${definitions[@]}"; do
    if [ $(grep -q "${tool^^}_VERSION=" "${tool}") ]; then
        version=$(grep "${tool^^}_VERSION=" "${tool}" | sed "s/.*=//")
    else
        version=$(grep "VERSION " "${tool}" | sed "s/.* //")
    fi;

    output_file="${tool}_${version}.sif";
    if [[ ! -f $output_file ]]; then
        echo "Building tool ${tool} to ${output_file}";
        singularity build --force "${output_file}" "${tool}";
        ln -sf "${output_file}" "${output_file%%_*}.sif"
    else
        echo "Tool ${tool} already exist, skipping...";
    fi;
done;
echo "Download pre built singularity containers";
for tool in ${!containers[@]}; do
    version=$(echo "${containers[$tool]}" | sed -r "s/.*://" | sed "s/--.*//")
    output_file="${tool}_${version}.sif";
    if [[ ! -f $output_file ]]; then
        echo "Downloading ${tool}";
        wget -O "${output_file}" "${containers[$tool]}"
        ln -s "${output_file}" "${output_file%%_*}.sif"
    else
        echo "Tool ${tool} already exist, skipping...";
    fi;
done

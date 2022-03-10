# JASEN
[![Docker Image CI](https://github.com/JD2112/JASEN/actions/workflows/docker-image.yml/badge.svg)](https://github.com/JD2112/JASEN/actions/workflows/docker-image.yml)
_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

## Setup
* `git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git`
* Edit `JASEN/nextflow.config`
* _`Optionally run: bash JASEN/container/safety_exports.sh USER PREFIX`_


## Singularity implementation
### Image creation
* Install Singularity (through conda or whatever)
* `cd JASEN/container && bash build_container.sh`

### Image execution
* `singularity exec -B JASEN_INSTALL_DIR:/external -B WORKDIR:/out IMAGE nextflow -C /external/nextflow.config run /JASEN/main.nf -profile local,singularity`


## Conda implementation
* Install Conda ( https://www.anaconda.com/distribution )
* Install nextFlow ( `curl -s https://get.nextflow.io | bash` )
* `bash JASEN/setup.sh`
* `nextflow run JASEN/main.nf -profile -local,conda`

# nextflow pipeline for typing and marker detection of bacteria

## Purpose

The pipeline is aimed at producing data useful for epidemiological and surveillance purposes. 
In v1 the pipeline is only tested using MRSA, but it should work well with
any bacteria having a good cgMLST scheme.

## Installation

Clone the pipeline repository with [nextflow-modules](https://github.com/Clinical-Genomics-Lund/nextflow-modules) submodule.

``` bash
git clone --recursive git@github.com:Clinical-Genomics-Lund/nextflow-modules.git
```

Install the database components required by the pipeline.

## How to use

Input files are defined in a csv file with the following format. All samples need to be of the same "type", meaning that they can be analyzed with the same analysis profile, defined in the nextflow config.

``` csv
id,read1,read2
p1,ALL504A259_122-78386_S1_R1_001.fastq.gz,ALL504A259_122-78386_S1_R2_001.fastq.gz
p2,ALL504A260_122-78386_S2_R1_001.fastq.gz,ALL504A260_122-78386_S2_R2_001.fastq.gz
p3,ALL504A261_122-78386_S3_R1_001.fastq.gz,ALL504A261_122-78386_S3_R2_001.fastq.gz
p4,ALL504A262_122-78386_S4_R1_001.fastq.gz,ALL504A262_122-78386_S4_R2_001.fastq.gz
p5,ALL504A263_122-78386_S5_R1_001.fastq.gz,ALL504A263_122-78386_S5_R2_001.fastq.gz
```

Start a new analsis with samples defined in `test.csv` using the staphylococcus_aureus profile.

``` bash
nextflow run -entry bacterial_default -profile staphylococcus_aureus -config configs/nextflow.trannel.config --csv=test.csv
```

## Components

### QC

Species detection is performed using [Kraken2](https://ccb.jhu.edu/software/kraken2/) together with [Bracken](https://ccb.jhu.edu/software/bracken/). 
The database used is a standard Kraken database built with 

``` bash
kraken2-build --standard --db $DBNAME
```

Low levels of Intra-species contamination or erronous mapping is removed using bwa and filtering away 
the heterozygous mapped bases. 

Genome coverage is estimated by mapping with [bwa mem](https://github.com/lh3/bwa) and using a bed file containing the cgMLST loci.

A value on the evenness of coverage is calculated as an [interquartile range](https://en.wikipedia.org/wiki/Interquartile_range).

### Epidemiological typing

For de novo asspembly [SPAdes](http://cab.spbu.ru/software/spades/) is used. [QUAST](http://cab.spbu.ru/software/quast/) 
is used for extraxting QC data from the assembly.

The cgMLST reference scheme used, is branched off [cgmlst.net](https://www.cgmlst.org/ncs/schema/141106/) 
At the moment this fork is not synced back with new allele numbers. For extracting alleles [chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki) 
is used. Number of missing loci is calculated and used as a QC parameter.

Traditional 7-locus MLST is calculated using [mlst](https://github.com/tseemann/mlst).

### Virulence and resistance markers

[ARIBA](https://github.com/sanger-pathogens/ariba) is used as the tool to detect genetic markes. 
The database for virulence markes is [VFDB](http://www.mgc.ac.cn/VFs/).

## Report and visualisation

The QC data is aggregated in a web service CDM (repo coming) and the cgMLST is visualized using a web service 
cgviz that is combined with [graptetree](https://github.com/achtman-lab/GrapeTree) for manipulating trees (repo coming).

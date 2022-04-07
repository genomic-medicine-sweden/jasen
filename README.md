<p align="center">
  <a href="https://github.com/genomic-medicine-sweden/JASEN">
    <img src="artwork/logo.png"/>
  </a>
</p>

_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

JASEN produces results for epidemiological and surveillance purposes. 
JASEN has been tested using MRSA, but should work well with any bacteria with a stable cgMLST scheme.

## Requirements

* Singularity
* Nextflow
* (Optional): Conda

## Conda deployment (self-contained)
* `git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git`
* `cd JASEN`
* `curl -s https://get.nextflow.io | bash`
* `bash deploy/deploy_conda/setup.sh`
* `bash deploy/deploy_references.sh` 
* `singularity remote login`
* `bash container/build_container.sh`

## Usage

### Simple self-test
``` bash
./nextflow run main.nf -entry bacterial_default -profile staphylococcus_aureus -config configs/nextflow.base.config --csv=assets/test_data/samplelist.csv
```

Start a new analysis with samples defined in `assets/test_data/samplelist.csv` using the staphylococcus_aureus profile.

Input files are defined in a csv file with the following format. All samples need to be of the same "type", meaning that they can be analyzed with the same analysis profile, defined in the nextflow config.

``` csv
id,read1,read2
p1,ALL504A259_122-78386_S1_R1_001.fastq.gz,ALL504A259_122-78386_S1_R2_001.fastq.gz
p2,ALL504A260_122-78386_S2_R1_001.fastq.gz,ALL504A260_122-78386_S2_R2_001.fastq.gz
p3,ALL504A261_122-78386_S3_R1_001.fastq.gz,ALL504A261_122-78386_S3_R2_001.fastq.gz
p4,ALL504A262_122-78386_S4_R1_001.fastq.gz,ALL504A262_122-78386_S4_R2_001.fastq.gz
p5,ALL504A263_122-78386_S5_R1_001.fastq.gz,ALL504A263_122-78386_S5_R2_001.fastq.gz
```

## Component Breakdown

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

# nextflow pipeline for typing and marker detection of bacteria

## Purpose

The pipeline is aimed at producing data useful for epidemiological and surveillance purposes. 
In v1 the pipeline is only tested using MRSA, but it should work well with
any bacteria having a good cgMLST scheme.

## Components
### QC

Species detection is performed using [Kraken2](https://ccb.jhu.edu/software/kraken2/) together with [Bracken](https://ccb.jhu.edu/software/bracken/). 
The database used is a standard Kraken database built with ```kraken2-build --standard --db $DBNAME```

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

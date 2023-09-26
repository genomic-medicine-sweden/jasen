<p align="center">
  <a href="https://github.com/genomic-medicine-sweden/JASEN">
    <img src="artwork/logo.png"/>
  </a>
</p>

_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

JASEN produces results for epidemiological and surveillance purposes.
JASEN has been developed for a small set of microbiota (primarily MRSA), but will likely work with any bacteria with a stable cgMLST scheme.

## Requirements

* Singularity
* Nextflow (`curl -s https://get.nextflow.io | bash`)

### Recommended

* Conda
* Singularity Remote Login

## Development deployment (self-contained)

### Copy code locally

```
git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git && cd JASEN
```

### Access to OCI registries (Optional)

```
singularity remote login
```

### Create singularity images. 
Running `build_container.sh` requires sudo priviledges.
```
cd container && sudo bash -i build_container.sh && cd ..
```

### Download references and databases using singularity. 
```
bash -i deploy/deploy_references_singularity.sh
```
Any errors produced during `deploy_references_singularity.sh` will hinder pipeline execution in unexpected ways.

## Configuration and test data

### Config 
Source: `configs/nextflow.base.config`

* Edit the `root` parameter in `configs/nextflow.base.config`
* Edit the `krakenDb`, `workDir` and `outdir` parameters in `configs/nextflow.base.config`
* Edit the `runOptions` in `configs/nextflow.base.config` in order to mount directories to your run

### Test data
Source: `assets/test_data/samplelist.csv`

* Edit the read1 and read2 columns in `assets/test_data/samplelist.csv`

## Setting up temp directories
Source: `~/.bashrc`

```
export SINGULARITY_TMPDIR="/tmp" #or equivalent
```

## Fetching databases

### Choose database
Choose between Kraken DB (64GB [Highly recommended]) or MiniKraken DB (8GB).
Or customize [your own](https://benlangmead.github.io/aws-indexes/k2).

### Download Kraken database

```
wget -O /path/to/kraken_db/krakenstd.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenstd.tar.gz
```

### Download MiniKraken database

```
wget -O /path/to/kraken_db/krakenmini.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenmini.tar.gz
```

## Updating databases

### Download Kraken database

```
bash /path/to/JASEN/assets/mlst_db/update_mlst_db.sh
```

## Usage

### Simple self-test

```
nextflow run main.nf -profile staphylococcus_aureus -config configs/nextflow.base.config --csv assets/test_data/samplelist.csv
```

### Usage arguments

| Argument type | Options                                | Required |
| ------------- | -------------------------------------- | -------- |
| -profile      | staphylococcus_aureus/escherichia_coli | True     |
| -entry        | bacterial_default                      | True     |
| -config       | nextflow.base.config                   | True     |
| -resume       | NA                                     | False    |
| --output      | user specified                         | False    |

### Input file format 

```csv
id,platform,read1,read2
p1,illumina,assets/test_data/sequencing_data/saureus_10k/saureus_large_R1_001.fastq.gz,assets/test_data/sequencing_data/saureus_10k/saureus_large_R2_001.fastq.gz
```

## Component Breakdown

### QC

* [Kraken2](https://ccb.jhu.edu/software/kraken2/): Species detection.
* [Bracken](https://ccb.jhu.edu/software/bracken/): Combined with Kraken2 for species detection.
* [bwa mem](https://github.com/lh3/bwa): Maps reads to cgMLST loci (demarcated by bed file) in order to estimate genome coverage. Low levels of Intra-species contamination or erroneous mapping is removed using bwa and filtering away the heterozygous mapped bases.
* [interquartile range](https://en.wikipedia.org/wiki/Interquartile_range): Calculates evenness of coverage.

### Assembly

* [SPAdes](http://cab.spbu.ru/software/spades/): De novo assembly for Ion Torrent.
* [SKESA](https://www.ridom.de/seqsphere/ug/v60/SKESA_Assembler.html): De novo assembly for Illumina.
* [QUAST](http://cab.spbu.ru/software/quast/): Extracts QC data (De novo assembly parameters) from the assembly.

### Epidemiological typing

* [chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki): Calculates cgMLST of extracted alleles decided by schema. Number of missing loci is calculated and used as a QC parameter.
* [cgmlst.net](https://www.cgmlst.org/ncs/schema/141106/): The cgMLST reference schema.
* [mlst](https://github.com/tseemann/mlst): Caculates traditional 7-locus MLST.

#### Supported profiles:

* `staphylococcus_aureus`
* `escherichia_coli`

#### Future profiles that will be supported:

* `klebsiella_pneumoniae`
* `mycobacterium_tuberculosis`

### Virulence and resistance markers

* [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/): Detects antimicrobial resistance genes as well as environmental and chemical resistance genes.
* [pointfinder](https://bitbucket.org/genomicepidemiology/pointfinder/src/master/): Combines with resfinder to detect variants.
* [virulencefinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/): Detects virulence genes.
* [amrfinderplus](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus): Detects antimicrobial resistance genes as well as environmental, chemical resistance and virulence genes.
* [resfinder_db](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/): Resfinder database.
* [pointfinder_db](https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/): Pointfinder database.
* [virulencefinder_db](https://bitbucket.org/genomicepidemiology/virulencefinder_db/src/master/): Virulencefinder database.

### Relatedness

* [sourmash](https://github.com/sourmash-bio/sourmash): Determine relatedness between samples.

## Report and visualisation

* [Bonsai](https://github.com/Clinical-Genomics-Lund/cgviz): Visualises JASEN outputs.
* [graptetree](https://github.com/achtman-lab/GrapeTree): Visualise phylogenetic relationship using cgmlst data.

## Tips

* Always run the latest versions of the bioinformatical software.
* Verify you have execution permission for JASENs `*.sif` images.
* Old Singularity versions may sporadically produce the error `FATAL: could not open image JASEN/container/*.sif: image format not recognized!`


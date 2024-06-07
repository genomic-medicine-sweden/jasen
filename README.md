<p align="center">
  <a href="https://github.com/genomic-medicine-sweden/jasen">
    <img src="artwork/logo.png"/>
  </a>
</p>

_Just Another System for Epityping NGS data_

>[!WARNING]
>**JASEN is in beta stage and the results are unverified. There is no guarantee that the pipeline can execute, output format consistency, or that it produces accurate results until there is an official 1.0 release.**

Jasen produces results for epidemiological and surveillance purposes.
Jasen has been developed for a small set of microbiota (primarily MRSA), but will likely work with any bacteria with a stable cgMLST scheme.

## Requirements

* [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-on-windows-or-mac)
* [JRE 8 - 21](https://www.java.com/en/download/manual.jsp)
* Nextflow (`curl -s https://get.nextflow.io | bash`)

### Recommended

* Conda
* Singularity Remote Login

## Usage

### Simple self-test

```
nextflow run main.nf -profile staphylococcus_aureus -config configs/nextflow.base.config --csv assets/test_data/samplelist.csv
```

#### Usage arguments

| Argument type | Options                                | Required |
| ------------- | -------------------------------------- | -------- |
| -profile      | **staphylococcus_aureus**, escherichia_coli, klebsiella_pneumoniae, mycobacterium_tuberculosis| True     |
| -config       | **configs/nextflow.base.config**, configs/nextflow.dev.config, configs/nextflow.hopper.config, configs/nextflow.ngp.config| True     |
| -entry        | bacterial_default                      | True     |
| --output      | User specified directory                         | False    |
| -resume       | Not applicable                                     | False    |


### Input file format 

```csv
id,platform,read1,read2
p1,illumina,assets/test_data/sequencing_data/saureus_10k/saureus_large_R1_001.fastq.gz,assets/test_data/sequencing_data/saureus_10k/saureus_large_R2_001.fastq.gz
```

### Update databases

#### Update MLST database

```
bash /path/to/jasen/assets/mlst_db/update_mlst_db.sh
```


## Installation

### Copy code locally

```
git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/jasen.git && cd jasen
```

### Create singularity images. 

The containers will be attempted to be built and downloaded as part of
the main Makefile (that is, when running `make install` in the main repo
folder).

```
cd container
make
```


### Download references and databases using singularity. 

First, make sure you stand in the `container` folder. Then run the `make` commands:

```
cd ..
make install
make check
```

Any errors produced during this step will hinder pipeline execution in
unexpected ways.

## Configuration

### Nextflow configuration
Source: `configs/nextflow.base.config`

* Edit the `root` parameter in `configs/nextflow.base.config`
* Edit the `krakenDb`, `workDir` and `outdir` parameters in `configs/nextflow.base.config`
* Edit the `runOptions` in `configs/nextflow.base.config` in order to mount directories to your run

When analysing Nanopore data:
* Edit the `ext.args` for Flye: specify genome size for the organism of interest with flag `--genome-size`
* Edit the `ext.seqmethod`for Flye depending on the input data
* Edit the `ext.args` for Medaka: specify the model with flag `-m`. Currently it is set to `r941_min_sup_g507`, but one should always set it based on how the data was produced. More about choosing the right model can be found [here](https://github.com/nanoporetech/medaka#models).

### Test data configuration
Source: `assets/test_data/samplelist.csv`

* Edit the read1 and read2 columns in `assets/test_data/samplelist.csv`

### Temporary directories configuration
Source: `~/.bashrc`

* Add the export line to `~/.bashrc`
* Change `SINGULARITY_TMPDIR` to `APPTAINER_TMPDIR` if you are using apptainer

```
export SINGULARITY_TMPDIR="/tmp" #or equivalent filepath to tmp dir
```

### Database configuration

#### Kraken database configuration
Choose between Kraken DB (64GB [Highly recommended]) or MiniKraken DB (8GB).
Or customize [your own](https://benlangmead.github.io/aws-indexes/k2).

##### Download standard Kraken database

```
wget -O /path/to/kraken_db/krakenstd.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenstd.tar.gz
```

##### (Alternatively) Download miniKraken database

```
wget -O /path/to/kraken_db/krakenmini.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenmini.tar.gz
```

#### Create TBProfiler database

##### Install jasentool

```
git clone git@github.com:ryanjameskennedy/jasentool.git && cd jasentool
pip install .
```

##### Create input csv that is used as tbdb input (composed of FoHM, WHO & tbdb variants)

```
jasentool converge --output_dir /path/to/jasen/assets/tbdb
```

##### Create tbdb (ensure tb-profiler is installed)

```
cd /path/to/jasen/assets/tbdb
tb-profiler create_db --prefix converged_who_fohm_tbdb
tb-profiler load_library converged_who_fohm_tbdb
```

##### Bgzip and index gms TBProfiler db

```
bgzip -c converged_who_fohm_tbdb.bed > /path/to/jasen/assets/tbprofiler_dbs/bed/converged_who_fohm_tbdb.bed.gz
tabix -p bed /path/to/jasen/assets/tbprofiler_dbs/bed/converged_who_fohm_tbdb.bed.gz
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
* [Flye](https://github.com/fenderglass/Flye/tree/flye): De novo assembly for Oxford Nanopore Technologies (ONT).
* [Medaka](https://github.com/nanoporetech/medaka): Creates consensus sequences from ONT data.

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

* [Bonsai](https://github.com/Clinical-Genomics-Lund/cgviz): Visualises jasen outputs.
* [graptetree](https://github.com/achtman-lab/GrapeTree): Visualise phylogenetic relationship using cgmlst data.

## Frequent issues / Tips

* Always run the latest versions of the bioinformatical software.
* Verify you have execution permission for jasens `*.sif` images.
* Old Singularity versions may sporadically produce the error `FATAL: could not open image jasen/container/*.sif: image format not recognized!`


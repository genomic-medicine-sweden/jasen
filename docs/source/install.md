# Installation

## Requirements

* Apptainer
* Nextflow (`curl -s https://get.nextflow.io | bash`)

**Recommended**

* Conda

## Development deployment (self-contained)

### Copy code locally

```bash
git clone --branch master              \\ 
    https://github.com/genomic-medicine-sweden/jasen.git && \\
cd jasen
```

### Installation requirements

**NOTE**: We assume that your OS has the following command-line tools installed in order for installation of JASEN:

```bash
unzip
gcc
zlib
```

### Create Apptainer images. 

The containers will be attempted to be built and downloaded as part of the main Makefile (that is, when running `make install` in the main repo folder).

```bash
cd containers && make
```

### Download references and databases using Apptainer. 

First, make sure your current working directory is in the main jasen folder (so if you cd:ed into the `container` folder before, you need to cd back to the main folder with `cd ..`). Then run the `install` make rule:

**NOTE**: Kraken and MLST databases need to be downloaded manually! Installation can be done independently for different species. Please see instructions below!

```bash
make install
```

Finally, run checks:

```bash
make check
```

Any errors produced during this step will hinder pipeline execution in unexpected ways.

### Species-specific installation

The following species are able be installed independently as to save time and disk usage:
 * saureus
 * ecoli
 * kpneumoniae
 * mtuberculosis

This is done by executing the following:

**NOTE**: `spyogenes` & `streptococcus` don't have any specific installation requirements, so `make update_databases` should suffice.

```bash
ORG="saureus"
make update_databases && make ${ORG}_all
```

## Configuration and test data

### Config 

Source: `nextflow.config`

* Edit the `root` parameter
* Edit the `workDir` and `outdir` parameters
* Edit the `use_kraken` parameter (default: false) and `kraken_db` to specify path to the database
* Edit the `use_hostile` parameter in `nextflow.config` in order to filter out human reads (default: false)
* Edit the `use_skesa` parameter (default: true) if you would like to use SPAdes instead of Skesa for assembly of short reads
* Edit the `target_sample_size` parameter in order to downsample reads
* Add  `runOptions` to apptainer/singularity profile in order to mount directories to your run, e.g. output folder, workdir (Example: `apptainer.runOptions = "--bind ${params.outdir} --bind ${params.workDir}"`)

When analysing Nanopore data:
* Edit the `ext.seqmethod` in `conf/modules.config` for Flye in case you are using older ONT data (default: --nano-hq, suitable for ONT data generated with R10 chemistry)

### Test data
Source: `assets/test_data/samplelist*.csv`

* For short reads produced with Illumina or IonTorrent technology, edit the `read1` and `read2` columns in `assets/test_data/samplelist.csv`
* For long reads produced with ONT technology, edit the `read1` column in `assets/test_data/samplelist_nanopore.csv`

## Setting up temp directories

Source: `~/.bashrc`

* Add the export line to `~/.bashrc`
* Change `SINGULARITY_TMPDIR` to `APPTAINER_TMPDIR` if you are using apptainer

```bash
export SINGULARITY_TMPDIR="/tmp" #or equivalent filepath to tmp dir
```

## Fetching/updating databases

**NOTE**: Both `kraken` and `mlst` require their databases to be downloaded **MANUALLY**

### Kraken

Choose between Kraken DB (64GB [Highly recommended]) or MiniKraken DB (8GB). Alternatively you can customize [your own](https://benlangmead.github.io/aws-indexes/k2).

#### Download Kraken database

```bash
wget -O /path/to/kraken_db/krakenstd.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenstd.tar.gz
```

#### Download MiniKraken database

```bash
wget -O /path/to/kraken_db/krakenmini.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenmini.tar.gz
```

### MLST databases (PubMLST & Blast)

**NOTE**: PubMLST DB requires users to have an account at [Bacterial Isolate Genome Sequence Database (BIGSdb)](https://pubmlst.org/bigsdb) in order to download the latest reported alleles. Here are the steps:
1. Register to all databases by clicking the `Database registrations`, check all, and register.
2. Create an API key under the `API keys` dropdown. 
3. Add them to `assets/mlstdb/update_mlstdb.sh` by editing `CLIENT_ID` and `CLIENT_SECRET`, or add the following to your `~/.bashrc`:
```
export PUBMLST_CLIENT_ID="<pubmlst_client_id>"
export PUBMLST_CLIENT_SECRET="<pubmlst_client_secret>"
```

#### Download/update MLST database

```bash
bash /path/to/jasen/assets/mlstdb/update_mlstdb.sh
```

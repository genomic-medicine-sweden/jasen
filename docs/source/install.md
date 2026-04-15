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
 * klebsiella
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

#### Batch mode (`kraken_batch`)

Set `use_kraken_batch = true` in `nextflow.config` to collect all samples into a single job and copy the Kraken database into `/dev/shm` (a RAM-backed filesystem) before classification. This loads the database once per pipeline run and eliminates disk I/O during classification, significantly reducing runtime for large sample batches.

**Requirements:** Each compute node running Kraken must have at least 80–100 GB of available `/dev/shm`. Check available space with:

```bash
df -h /dev/shm
```

If your nodes do not meet this requirement, leave `use_kraken_batch = false` (default) to run kraken2 per sample.

### cgMLST database (BIGSdb Pasteur setup)

**NOTE**: The *Klebsiella* cgMLST schema is hosted on [BIGSdb Pasteur](https://bigsdb.pasteur.fr/) and requires API credentials to download. Here are the steps:
1. Request an API key by following the instructions at [https://bigsdb.pasteur.fr/requesting-api-key/](https://bigsdb.pasteur.fr/requesting-api-key/).
2. Copy the client credentials template:
```
cp assets/.bigsdb_tokens/client_credentials.template assets/.bigsdb_tokens/client_credentials
```
3. Edit `client_id` and `client_secret` in the `assets/.bigsdb_tokens/client_credentials` file.
```
[Pasteur]
client_id = insert_pasteur_client_id
client_secret = client_id = insert_pasteur_client_secret
```
4. To download the raw cgMLST alleles from BIGSdb Pasteur, run:
**NOTE**: This target must be run manually and is **not** part of `make install`. It requires OAuth credentials to be configured as described above.
```bash
make klebsiella_download_cgmlst_schema
```
5. After downloading, re-reference the alleles by running:
```bash
make klebsiella_prep_cgmlst_schema
```

### MLST databases (PubMLST & BLAST)

**NOTE**: PubMLST DB requires users to have an account at [Bacterial Isolate Genome Sequence Database (BIGSdb)](https://pubmlst.org/bigsdb) in order to download the latest reported alleles. Here are the steps:
1. Register to all databases by clicking the `Database registrations`, check all, and register.
2. Create an API key under the `API keys` dropdown.
3. Add your credentials to your `~/.bashrc`:
```bash
export PUBMLST_CLIENT_ID="<pubmlst_client_id>"
export PUBMLST_CLIENT_SECRET="<pubmlst_client_secret>"
export PASTEUR_CLIENT_ID="<pasteur_client_id>" # From BIGSdb Pasteur setup
export PASTEUR_CLIENT_SECRET="<pasteur_client_secret>" # From BIGSdb Pasteur setup
```

#### Download/update MLST database per species

Run the token setup step first, then the database build step. Both steps require the `PUBMLST_CLIENT_ID` and `PUBMLST_CLIENT_SECRET` (PubMLST schemas) or `PASTEUR_CLIENT_ID` and `PASTEUR_CLIENT_SECRET` (Pasteur schemas) environment variables.

**S. aureus**
```bash
make setup_saureus_mlstdb_token
make update_saureus_mlstdb
```

**S. pyogenes**
```bash
make setup_spyogenes_mlstdb_token
make update_spyogenes_mlstdb
```

**E. coli achtman**
```bash
make setup_ecoli_achtman_mlstdb_token
make update_ecoli_achtman_mlstdb
```

**E. coli pasteur** (needs BIGSdb Pasteur setup)
```bash
make setup_ecoli_pasteur_mlstdb_token
make update_ecoli_pasteur_mlstdb
```

**Klebsiella** (needs BIGSdb Pasteur setup)
```bash
make setup_klebsiella_mlstdb_token
make update_klebsiella_mlstdb
```

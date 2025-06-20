# genomic-medicine-sweden/jasen: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Fixed

### Changed

## [1.1.1]

### Added

### Fixed

- Fixed chewbbaca arg bug

### Changed

- Changed computational resources for hostile

## [1.1.0]

### Added

- Added missing tools to S. aureus subworkflow tool list
- Added spyogenes bed file representing cgmlst targets
- Added `prodigal` image download to container `Makefile`
- Added module for mapping long reads with Minimap2
- Added Minimap2 refence index `referenceGenomeMmi` for all species
- Added `samtools_coverage` to get stats from mapping long reads to the reference and assembly
- Added `minimap2_align_assembly`
- Added `minimap2_index` for assembly indexing
- Added `when` operator to `samtools` for `nanopore`
- Added `assay` to yaml output
- Added `staphylococcus_nrl` profile
- Added `release_life_cycle` via profiles
- Added `create_prp_yaml.py` to bin
- Added `spatyper` and `sccmec` to yaml
- Added `gambitcore` as a module
- Added `apptainer` and `singularity` profiles

### Fixed

- Fixed stub-run for S.aureus ONT workflow
- Fixed typing error in Makefile
- Fixed checking out pipeline code with submodules
- Fixed genome downloading via `bin/download_ncbi.py` with timed retries
- Fixed dirty submodules in `.gitmodules`
- Fixed dubious ownership bug where finder dbs retrieve commit ID

### Changed

- Changed spyogenes genome from `GCF_000006785.2` to `GCF_005164585.1`
- Updated `check_taxon` method to include spyogenes
- Changed ptf downloading to generation of ptf via prodigal for `ecoli` & `spyogenes`
- Renamed index created by BWA from `referenceGenomeIdx` to `referenceGenomeFai`
- Renamed `bwa_mem_dedup` to `bwa_mem_assembly`
- Split `samtools_sort` into two modules
- Removed `abritamr` as it's not used
- Changed dirname `container` to `containers`
- Changed dirname `configs` to `conf`
- Changed dirname `nextflow-modules` to `modules`
- Changed importing of `spades`
- Updated `Makefile` re dir name changes
- Updated subworkflows
- Grouped most workflows into one (`bacterial_general.nf`) to remove duplicate code
- Changed to have only one main config `nexflow.config`
- Changed `cmd` module name to `cdm`
- Updated bonsai-prp to v1.3.1
- Changed indentation structure
- Moved `platform` to config via `params.platform` 
- Changed `hostile` io
- Updated docs regarding restructuring
- Changed `prp` sub commands
- Updated CI GA workflow re `container_dir`
- Removed `check-and-reinit-git-submodules` from CI GA workflow
- Removed `kma` submodule
- Updated finder submodules (`virulencefinder_db`, `resfinder_db`, `pointfinder_db`, `serotypefinder_db`)
- Updated `resfinder` version
- Changed memory allowance for hostile
- Updated docs regarding updating finder dbs

## [1.0.0]

### Added

- Light installation doc updates
- Updated NGP config file
- Added ska2 module process (`ska_build`)
- Added ska2 to `bacterial_base.nf`
- Added emmtyper module process (`emmtyper`)
- Added emmtyper to `Makefile`
- Added downloading of Streptococcus & Streptococcus pyogenes files to `makefile`
- Added `get_taxon` to methods
- Added `streptococcus` and `streptococcus_pyogenes` workflows and profiles to configs
- Added ska filepath to yaml
- Added optional read downsampling using seqtk
- Added `tbGradingRulesBed` to IGV track
- Added module `nanoplot` to check quality of raw reads from ONT
- Added module `fastq` to check qc of strep samples
- Added module `hostile` to remove human reads
- Added module `sccmec` for typing of SCCmec cassettes in assemblies of Staphylococcus species
- Added `when` operators to all modules that aren't in `bacterial_base.nf`
- Added additional arguments for Resfinder when analysing Nanopore data
- Added downloading of hostile index to `Makefile`
- Added sw docs for `kpneumonia`, `spyogenes`, & `streptococcus`
- Added module `spatyper` for typing of Spa gene in Staphylococcus aureus

### Fixed

- Interrupted installation now hard resets E.Coli. to avoid fragmented installation
- Installation target file for streptococcus changed, to resolve interrupted installation
- Included htslib image to natively support tabix and bgzip
- Fixed a bug that prevented pipeline from correctly guessing misspelled profiles
- Fixed tbprofiler related installation rules in `Makefile`
- All modules can be called in one workflow -> profile determines whether they are actually executed using `when` statement
- Empty channels fixed
- Chewbbaca collection of input fixed
- Fixed resfinder `--species` arg
- Fixed `nextflow.hopper.config` `symlinkDir`
- Removed serotypefinder from saureus workflow
- Fixed jasen running only on the first row/sample in csv
- Fixed channel problem by changing `Channel.of([])` to `Channel.value([])`
- Fixed medaka by changing `--threads` to `-t` in all the configs
- Fixed `nextflow.hopper.config` re singularity image path
- Fixed sccmec version file
- Fixed handling of ONT samples re fastqc & postalignqc
- Fixed io of spatyper and turned on postalignqc for ONT
- Fixed io of resfinder for all workflows but saureus
- Fixed input array for postalignqc

### Changed

- Removed sudo requirement from apptainer container creation (only required for older versions of apptainer) to streamline installation
- Updated tbdb submodule
- Moved taxon related methods to methods directory
- Changed spyogenes genome from GCF_900475035.1 to GCF_000006785.2
- Changed `containerDir` filepath for Lund configs
- Update `configs/nextflow.dev.config` root fpath
- Changed the freebayes output directory
- Remove `copy_to_cron` module
- Updated from Singularity v3.8.7 to Apptainer v1.3.6 in the CI pipeline
- Updated PRP to version 0.11.4
- Changed SerotypeFinder input from reads to assembly
- Changed variable formatting of modules
- Changed ska & sourmash filepath to symlink dir in `nextflow.hopper.config` & `nextflow.dev.config`
- Updated Kraken db filepath in `nextflow.hopper.config` & `nextflow.dev.config`
- Changed `staphylococcus_aureus_all` to `saureus_all` in `Makefile`
- Moved all `when` commands to configs
- Updated `fastqc` & `sccmec` mem settings
- Changed model that is used in `medaka_consensus` to bacterial model (using `--bacteria` argument) 

## [0.9.0]

### Added

- Added long-read test data (S. aureus)
- Added `samplelist_nanopore.csv` for running long-read test data
- Added location of documentation to `README` 
- Added `cdmDir` to config
- Added NanoPlot module
- Added process for adding IGV annotation tracks with PRP.
- Updated how `mycobacterium_tuberculosis` workflow adds IGV annotation tracks.

### Fixed

- Fixed `--qc` argument filepath to be full filepath to output
- Fixed tbprofiler url in container Makefile
- Fixed TB installation steps in main Makefile

### Changed

- Updated TbProfiler to version 6.3
- Updated PRP to version 0.10.0
- Removed delly annotation
- Updated vcf args in prp module

## [0.8.0]

### Added

- Added ShigaPass
- Added mlstBlastDb to mlst
- Added full path for bam and vcf filepaths 
- Added bam and bai to bonsai input for `staphylococcus_aureus`, `escherichia_coli` & `klebsiella_pneumoniae`
- Added `bamDir` and `vcfDir` to config params
- Added run `bwa_mem` from only when profile is not `mycobacterium_tuberculosis`
- Added `sample_name` to `get_seqrun_meta`
- Automatically publish the pipeline documentation to read the docs.

### Fixed

- ShigaPass URL fixed
- Fixed qc channel regarding `mycobacterium_tuberculosis`
- Fixed bwa output file bug and stub
- Fixed README re jasentool cmds
- Fixed `get_seqrun_meta` if statements
- Fixed getting some software versions
- Fixed tb workflow bug

### Changed

- Fixed output format for tbprofiler
- Removed `samtools_sort_ref` from configs
- Changed `--symlink_dir` arg for prp
- Changed sampleName to sampleID
- Updated `flye`, `freebayes`, `mask`, `medaka`, `post_align_qc`, `skesa` & `spades`  output filenames
- Updated bonsai-prp to v0.9.3
- Updated the pipeline documentation.
- Removed sudo from make (deprecated)
- Updated NGP config for new hardware
- Updated tbprofiler memory allocation

## [0.7.0]

### Added

- Added two modules for assembly of long-read Nanopore data: Flye and Medaka
- Added `referenceGenome` to `prp create-bonsai-input`
- created `.fai` files from all genomes in `Makefile`
- Added `samtools_sort_ref` to tb workflow with 4GB memory (may need more)
- Added more cpus to tbprofiler

### Fixed

- Added nanopore to all workflows
- Sort and index tbprofiler bam output

### Changed

- Updated saureus genome from `NC_002951.2` to `GCF_000012045.1`
- Updated ecoli genome from `NC_000913.3` to `GCF_000005845.2`
- Updated mtuberculosis genome from `NC_000962.3` to `GCF_000195955.2`
- Updated kpneumoniae genome from `NC_016845.1` to `GCF_000240185.1`
- Updated tbprofiler to v6.2.0
- Updated `download_ncbi.py` to include `.gff` files
- Updated bonsai-prp to v0.8.3
- Removed FoHM variant duplicates

## [0.6.0]

### Added

- Added converged_who_fohm_tbdb.csv
- Added guide to create tbdb
- Added `sequencing_run` and `lims_id` to output
- Added `devMode` flag
- Added `lims_id` & `sequencing_run` to `meta` module
- Added `annotate_delly` module
- Added prp versions
- Added how-to guide on how to create bed file, bgzipped format and index for `annotate_delly` input
- Prepare general `Makefile` for incorporation of tbprofiler v6.1.0 to auto create above files and create new tbdb
- Add `create_yaml` module for upload to bonsai

### Fixed

- Args for cdm

### Changed

- Run with converged/merged db
- Publish bam & bai from tbprofiler
- cronCopy set true in hopper config
- bonsai-prp version upgraded from v0.5.0 to v0.6.0
- Renamed `nextflow.hopper.config` to `nextflow.dev.config` for hopper development
- bonsai-prp version upgraded from v0.6.0 to v0.7.0
- bonsai-prp version upgraded from v0.7.0 to v0.7.1
- Changed blastDb to pubMlstDb re mlst

## [0.5.0]

### Added

- Added a GitHub workflow to run a basic CI pipeline.
- Build prp as singularity image from dockerhub in Makefile
- Chewbacca and virulencefinder sif fetched from galaxyproject
- Added species to amrfinderplus
- Added getSpeciesTaxonName to amrfinderplus
- Add serotypefinder to Makefiles, workflows, modules & configs

### Fixed

- Fixed config not pointing to the new lowercase repo name: jasen (instead of JASEN)

### Changed

- Changed to use full file paths in include statements for better navigation in text editors.
- Upgraded the Skesa (container) to v2.5.1 to fix ownership issue with /tmp folder
- Changed pythonScripts.sif filename to bonsai-prp.sif
- Postalignqc added to prp and move to prp module
- Makefile converts bonsai-prp from docker image to singularity
- Makefile pulls all containers from online respositories
- Configs can also pull images from galaxy project or dockerhub
- Move post_align_qc to prp module

## [0.4.0]

### Added

- Added [tbdb](https://github.com/jodyphelan/tbdb) as a submodule.
- Downloaded new catalogue from WHO and convert to csv required for tbdb creation using [jasentool](https://github.com/ryanjameskennedy/jasentool).
- DB creation included using `Makefile`.
- Updated configs.
- Updated `workflows/mycobacterium_tuberculosis.nf`.
- Update TBProfiler version to v5.0.1.
- Publish delly output when running tbprofiler.
- Remove Tbprofiler run using WHO database until TBProfiler issue #313 resolved.

### Fixed


### Changed

- Updated PRP to verion 0.3.0

## [0.3.0](https://github.com/genomic-medicine-sweden/JASEN/tag/v0.3.0)

### Added

- _E.coli_ STX typing using Virulencefinder

### Changed

- Moved PRP to seperate repo
- Updated PRP to v0.2.0

## [0.2.0](https://github.com/genomic-medicine-sweden/JASEN/tag/0.2.0-rc)

### Added

- The pipeline now supports IonTorrent input data
- Both Single-End and Paired-End input is now supported
- Test suite is entirely containerised within Singularity
- Versions of all used software is conveniently written into a single file
- Assemblies are automatically cleaned up
- AMRFinder has replaced ARIBA
- SKESA de novo assembly added
- Readme contains simple usage instructions and more detailed breakdown of components
- External nextflow-modules are now fully integrated

### Changed

- Extensive Nextflow code refactoring
- Memory and CPU usage defaults optimised for all modules
- Optimisation of module invocations to reduce overhead
- PRP identifiers updated
- Updated all utilised databases to latest version
- Various publishDir updates

### Fixed

- Resolved multiple BLAST database referencing errors
- Singularity images are now umasked properly
- Minimum recommended RAM expanded to avoid runtime crashes
- Fixed 'failed to publish file' error
- MLST PRP input issues fixed
- Metadata is now properly saved

### Removed

- Removed legacy modules (notably ARIBA)
- Reduction to only two software requirements; Singularity and Nextflow

## [0.1.0](https://github.com/genomic-medicine-sweden/JASEN/tag/0.1.0-beta)

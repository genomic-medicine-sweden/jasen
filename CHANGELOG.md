# genomic-medicine-sweden/jasen: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added a GitHub workflow to run a basic CI pipeline.
 - Build prp as singularity image from dockerhub in Makefile
 - Chewbacca and virulencefinder sif fetched from galaxyproject
 - Added species to amrfinderplus
 - Added getSpeciesTaxonName to amrfinderplus

### Fixed

 - Fixed config not pointing to the new lowercase repo name: jasen (instead of JASEN)

### Changed

 - Changed to use full file paths in include statements for better navigation in text editors.
 - Upgraded the Skesa (container) to v2.5.1 to fix ownership issue with /tmp folder
 - Changed pythonScripts.sif filename to bonsai_prp.sif
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

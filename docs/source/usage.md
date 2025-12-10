# Usage

## Simple self-test

with Illumina test data:
```bash
nextflow run main.nf                                      \
        -profile staphylococcus_aureus,illumina,apptainer \
        -config nextflow.config                           \
        --csv assets/test_data/samplelist.csv
```

with ONT test data:
```bash
nextflow run main.nf                                      \
        -profile staphylococcus_aureus,nanopore,apptainer \
        -config nextflow.config                           \
        --csv assets/test_data/samplelist_nanopore.csv
```
Remember to edit the paths to the test file(s) in the samplelist.

## Usage arguments

| Argument type       | Options                                                                                                | Required |
| ------------------- | ------------------------------------------------------------------------------------------------------ | -------- |
| -profile (species)  | staphylococcus_aureus/escherichia_coli/mycobacterium_tuberculosis/klebsiella_pneumoniae/streptococcus_pyogenes/streptococcus | True     |
| -profile (platform) | illumina/nanopore/iontorrent                                                                           | True     |
| -profile (RLS)      | development/diagnostic/validation                                                                      | False    |
| -config             | nextflow.config                                                                                        | True     |
| -resume             | /NA                                                                                                    | False    |
| --output            | user-specified                                                                                         | False    |

RLS = Release life cycle (default: diagnostic)

## Input file format 

```{csv-table} Example of a *samplelist* input file in CSV format.
:header-rows: 1

id,platform,read1,read2
p1,illumina,assets/test_data/sequencing_data/saureus_10k/saureus_large_R1_001.fastq.gz,assets/test_data/sequencing_data/saureus_10k/saureus_large_R2_001.fastq.gz
```

```{csv-table} Example of a *samplelist* input file in CSV format.
:header-rows: 1

id,platform,read1
test_mecA,nanopore,assets/test_data/sequencing_data/saureus_1k/saureus_1k_ont.fastq.gz
```

## Downsampling reads

There are an option to use [seqtk](https://github.com/lh3/seqtk) downsample the number of for a sample as a preprocessing step before all other analyses. This can be useful if a sample was sequenced too deeply, as extreme sequencing depth can causes issues with *de-novo* assemblies.

Activate downsampling by setting the parameter `target_sample_size` to the either the desired number of reads or the fraction of reads to include in the config.

## Removing Human reads

There are an option to use [hostile](https://github.com/bede/hostile) to filter human reads from further analyses. This can be useful if a sample has been contaminated, which could cause issues with *de-novo* assemblies.

Activate human read depletion by setting the parameter `use_hostile` to `true` in the config.

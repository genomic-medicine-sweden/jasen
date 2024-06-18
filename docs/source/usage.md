# Usage

## Simple self-test

```bash
nextflow run main.nf                         \
        -profile staphylococcus_aureus       \
        -config configs/nextflow.base.config \
        --csv assets/test_data/samplelist.csv
```

## Usage arguments

| Argument type | Options                                | Required |
| ------------- | -------------------------------------- | -------- |
| -profile      | staphylococcus_aureus/escherichia_coli | True     |
| -entry        | bacterial_default                      | True     |
| -config       | nextflow.base.config                   | True     |
| -resume       | NA                                     | False    |
| --output      | user specified                         | False    |

## Input file format 

```{csv-table} Example of a *samplelist* input file in CSV format.
:header-rows: 1

id,platform,read1,read2
p1,illumina,assets/test_data/sequencing_data/saureus_10k/saureus_large_R1_001.fastq.gz,assets/test_data/sequencing_data/saureus_10k/saureus_large_R2_001.fastq.gz
```
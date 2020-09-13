#!/usr/bin/env nextflow

params.in="$baseDir/myraw.bcf"
snpcalling_output=file(params.in)


process bcf_vcf_tsv {

  input:
  file bcf_file from snpcalling_output

  output:
  file 'vcftools.recode.vcf' into bcf_to_vcf_output

  script:
  """
  bcftools view ${bcf_file} > vcftools.recode.vcf 
  """

}


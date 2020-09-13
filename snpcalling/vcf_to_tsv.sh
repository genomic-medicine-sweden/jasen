#!/bin/bash 

#vk vcf2tsv wide --print-header myraw.vcf 

gatk VariantsToTable -V myraw.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F I16 -F QS -F MQ0F -GF PL -O output.tsv


#!/bin/bash

af_fname=AFs.csv
bam_fname=$1

# write relevant positions to a dummy VCF
cp vcf_header.txt vars.vcf
sed '1d' "$af_fname" | tr '_' ',' | awk -F',' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $1, ".", $2, $3, ".", ".", "AF=1"}' |
    sort -k2,2 -n >> vars.vcf

# now run freebayes and format the output
freebayes -f refgenome.fa "$bam_fname" \
    --variant-input vars.vcf \
    --only-use-input-alleles |
    bcftools norm -f refgenome.fa -m - |
    bcftools query -f '%POS\_%REF\_%ALT\t[%GT\t%DP]\n'

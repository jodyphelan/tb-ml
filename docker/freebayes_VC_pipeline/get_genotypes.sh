#!/bin/bash

af_file=AFs.csv
bam_file=aligned_reads

# write relevant positions to a dummy VCF and a dummy BED file
cp vcf_header.txt vars.vcf
sed '1d' "$af_file" | tr '_' ',' | awk -F',' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $1, ".", $2, $3, ".", ".", "AF=1"}' |
    sort -k2,2 -n >> vars.vcf
sed '1d' "$af_file" | tr '_' ',' | awk -F',' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $1-1, $1, ".", ".", $2, $3, ".", "AF=1"}' |
    sort -k2,2 -n > vars.bed

# extract the reads overlapping with the variants of interest
# (index the BAM first)

samtools index "$bam_file"
time samtools view -bML vars.bed "$bam_file" > extracted.bam

# # now run freebayes and format the output
# samtools faidx refgenome.fa
# freebayes -f refgenome.fa extracted.bam \
#     --variant-input vars.vcf \
#     --only-use-input-alleles |
#     bcftools norm -f refgenome.fa -m - |
#     bcftools query -f '%POS\_%REF\_%ALT\t[%GT\t%DP]\n'

/bin/bash
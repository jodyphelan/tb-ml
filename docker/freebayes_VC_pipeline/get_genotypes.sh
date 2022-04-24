#!/bin/bash

# make sure that we are in `/data` (in case the container is run via singularity)
cd /data
bam_fname=aligned_reads
# use a default DP threshold of 10
dp_threshold=${1:-10}

# get relevant positions from STDIN and write them to a dummy VCF and a dummy BED file
mv vcf_header.txt vars.vcf
tee target_vars_AF.csv | sed '1d' | awk -F',' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $1, ".", $2, $3, ".", ".", "AF=1"}' |
    sort -k2,2 -n >> vars.vcf
grep -v "#" vars.vcf | awk -F'\t' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $2-1, $2, ".", ".", $3, $4, ".", "AF=1"}' |
    sort -k2,2 -n > vars.bed

# extract the reads overlapping with the variants of interest
# (we have to index the BAM first)
samtools index "$bam_fname"
samtools view -bML vars.bed "$bam_fname" -T refgenome.fa >extracted.bam

(
    # print the header for the result
    echo 'POS,REF,ALT,GT,DP'
    # now run freebayes and format the output
    freebayes -f refgenome.fa extracted.bam \
        --variant-input vars.vcf \
        --only-use-input-alleles |
        bcftools norm -f refgenome.fa -m - 2> >(grep -v ^Lines) |
        bcftools query -f '%POS,%REF,%ALT,[%GT,%DP]\n'
) | python process_variants.py "$dp_threshold"

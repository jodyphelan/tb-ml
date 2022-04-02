#!/bin/bash
# cd "$HOME"

af_fname=$1
bam_fname=$2
# if no DP threshold is given, use 10 as default value
DP_threshold="${3:-10}"

# write relevant positions to a dummy VCF and a dummy BED file
cp vcf_header.txt vars.vcf
sed '1d' /data/"$af_fname" | tr '_' ',' | awk -F',' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $1, ".", $2, $3, ".", ".", "AF=1"}' |
    sort -k2,2 -n >>vars.vcf
sed '1d' /data/"$af_fname" | tr '_' ',' | awk -F',' \
    'BEGIN{OFS="\t"}
    {print "Chromosome", $1-1, $1, ".", ".", $2, $3, ".", "AF=1"}' |
    sort -k2,2 -n >vars.bed

# extract the reads overlapping with the variants of interest
# first, index the BAM (create symlink to make things easier)
ln -s /data/"$bam_fname" .
samtools index "$bam_fname"
samtools view -bML vars.bed "$bam_fname" >tmp
mv tmp "$bam_fname"

# now run freebayes and format the output
freebayes -f refgenome.fa "$bam_fname" \
    --variant-input vars.vcf \
    --only-use-input-alleles |
    bcftools norm -f refgenome.fa -m - |
    bcftools query -f '%POS\_%REF\_%ALT\t[%GT\t%DP]\n' >new_vars.tsv

python sanitise_variants.py "$af_fname" "$DP_threshold" 

/bin/bash

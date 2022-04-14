FROM continuumio/miniconda3:latest

RUN conda install -c bioconda -c conda-forge -y \
    pandas=1.4.1 \
    freebayes=1.3.6 \
    bcftools=1.15 \
    samtools=1.14

RUN mkdir /data
WORKDIR /data

# copy refgenome header for dummy VCF
COPY MTB-h37rv_asm19595v2-eg18.fa ./refgenome.fa
COPY vcf_header.txt .

# copy entrypoint script and the python helper script
COPY get_genotypes.sh /
COPY process_variants.py /

ENTRYPOINT [ "/bin/bash", "/get_genotypes.sh" ]
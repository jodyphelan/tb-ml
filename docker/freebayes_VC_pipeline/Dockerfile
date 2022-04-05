FROM continuumio/miniconda3:latest

RUN conda install -c bioconda -c conda-forge -y \
    freebayes=1.3.6 \
    bcftools=1.15

RUN mkdir /data
WORKDIR /data

# copy refgenome header for dummy VCF
COPY MTB-h37rv_asm19595v2-eg18.fa ./refgenome.fa
COPY vcf_header.txt .

# copy entrypoint script
COPY get_genotypes.sh /

ENTRYPOINT [ "/bin/bash", "/get_genotypes.sh" ]
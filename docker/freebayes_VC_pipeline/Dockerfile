# FROM biocontainers/biocontainers:latest
FROM mambaorg/micromamba:0.22.0

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root
RUN mkdir /data
WORKDIR /data

# copy refgenome header for dummy VCF
COPY MTB-h37rv_asm19595v2-eg18.fa ./refgenome.fa
COPY vcf_header.txt .

# copy entrypoint script
COPY get_genotypes.sh .

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "/bin/bash", "get_genotypes.sh" ]
# ENTRYPOINT [ "/bin/bash" ]
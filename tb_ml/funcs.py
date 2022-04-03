from glob import glob
from uuid import uuid4
import pandas as pd
import numpy as np
import subprocess
import argparse
import os

from pyrsistent import b


def get_cli_args() -> argparse.Namespace:
    pass


def get_positions(AFs: pd.Series) -> list[int]:
    """
    Extracts the genomic positions of interest from the AF series; assumes
    the index entries to have the format: 'POS_REF_ALT'.
    """
    pos = [int(x.split("_")[0]) for x in AFs.index]
    return pos


def get_genotypes(
    bam_file: str,
    ref_file: str,
    AFs: pd.Series,
    subset_bam: bool = True,
    DP_threshold: int = 10,
) -> pd.Series:
    """
    1. Writes a dummy vcf with the variants
    2. Subsets the bam file to only include overlapping reads
    3. Uses freebayes to produce calls for those specific positions.
    """
    prefix = str(uuid4())
    tmp_vcf_file = f"{prefix}.vcf"
    tmp_bed_file = f"{prefix}.bed"
    tmp_bam_file = f"{prefix}.bam"
    # write the dummy VCF
    with open(tmp_vcf_file, "w") as outfile:
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write(
            "##reference=/home/jody/refgenome/" "MTB-h37rv_asm19595v2-eg18.fa\n"
        )
        outfile.write("##contig=<ID=Chromosome,length=4411532>\n")
        outfile.write(
            "##INFO=<ID=AF,Number=A,Type=Float,Description="
            '"Estimated allele frequency in the range (0,1]">\n'
        )
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for var in AFs.index:
            pos, ref, alt = var.split("_")
            outfile.write(f"Chromosome\t{pos}\t.\t{ref}\t{alt}\t.\t.\tAF=1\n")
    if subset_bam:
        # extract aligned reads overlapping with the variants of interest
        subprocess.run(f"vcf2bed < {tmp_vcf_file} > {tmp_bed_file}", shell=True)
        subprocess.run(
            f"samtools view -bML {tmp_bed_file} {bam_file} > {tmp_bam_file}", shell=True
        )
        bam_file = tmp_bam_file
    # variant calling with freebayes and output formateed with bcftools
    VC_result = subprocess.run(
        f"freebayes -f {ref_file} {bam_file} --variant-input {tmp_vcf_file} "
        "--only-use-input-alleles  "
        f"| bcftools norm -f {ref_file} -m - "
        r"| bcftools query -f '%POS\_%REF\_%ALT\t[%GT\t%DP]\n'",
        shell=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    # extract the variants and bring into a suitable form
    variants = pd.Series(VC_result.split("\n")).str.split("\t", expand=True)
    variants.columns = ["varID", "GT", "DP"]
    variants = variants.set_index("varID")
    variants["GT"] = variants["GT"].apply(lambda x: x[0])
    # declare variants with DP < threshold as non-calls and replace all
    # noncalls with the corresponding AF values
    noncalls_idx = [
        i
        for i, row in variants.iterrows()
        if row["GT"] == "." or row["DP"] == "." or int(row["DP"]) < DP_threshold
    ]
    variants.loc[noncalls_idx, "GT"] = AFs[noncalls_idx]
    # add the allele frequencies for variants not found in the variant
    # calling pipeline to ensure matching dimensions for the prediction model
    variants = pd.concat(
        (variants["GT"], AFs[[x for x in AFs.index if x not in variants.index]])
    )
    # make sure the order is as expected by the model
    variants = variants[AFs.index]
    # clean up temporary files
    for f in glob(f"{prefix}*"):
        os.remove(f)
    return variants


def load_container(container_tar_path: str) -> str:
    # load the container and return name of the generated image
    p = subprocess.run(
        ["docker", "load", "-i", container_tar_path], capture_output=True, text=True
    )
    if p.returncode:
        print(f"ERROR loading docker image '{container_tar_path}':")
        print(p.stderr)
        exit(1)
    img_name = p.stdout.split("Loaded image: ")[1].strip()
    return img_name


def run_VC_container(
    VC_container_tar_path: str, bam_file: str, af_file: str
) -> pd.Series:
    # load the container
    img_name = load_container(VC_container_tar_path)
    # run the container
    af_path = af_file if os.path.isabs(af_file) else f'{os.getcwd()}/{af_file}'    
    bam_path = bam_file if os.path.isabs(bam_file) else f'{os.getcwd()}/{bam_file}'    
    p = subprocess.run(
        [
            "docker",
            "run",
            "--mount",
            f"type=bind,source={af_file},target=/data/AFs.csv",
            "--mount",
            f"type=bind,source={bam_file},target=/data/aligned_reads",
            img_name,
        ],
        capture_output=True,
        text=True,
    )
    if p.returncode:
        print(f"ERROR running variant-calling pipeline with '{VC_container_tar_path}':")
        print(p.stderr)
        exit(1)
    # reformat the result
    variants = pd.Series(p.stdout.strip().split("\n")).str.split("\t", expand=True)
    variants.columns = ["varID", "GT", "DP"]
    variants = variants.set_index("varID")
    variants["GT"] = variants["GT"].apply(lambda x: x[0])
    variants = variants.replace(".", np.nan).astype(float)
    return variants


def sanitise_variants(
    variants: pd.Series, AFs: pd.Series, DP_threshold: int = 10
) -> pd.Series:
    new_vars = variants.copy()
    # declare variants with DP < threshold as non-calls and replace all
    # noncalls with the corresponding AF values
    new_vars.loc[(new_vars["DP"] < DP_threshold) | new_vars["DP"].isna(), "GT"] = np.nan
    new_vars["GT"].fillna(AFs, inplace=True)
    # add the allele frequencies for variants not found in the variant
    # calling pipeline to ensure matching dimensions for the prediction model
    new_vars = pd.concat((new_vars["GT"], AFs[AFs.index.difference(new_vars.index)]))
    # make sure the order is as expected by the model
    new_vars = new_vars[AFs.index]
    return new_vars


def run_prediction_container(pred_container_tar_path: str, vars: pd.Series) -> bool:
    # load the container
    img_name = load_container(pred_container_tar_path)
    # run the container to get the prediction
    p = subprocess.run(
        ["docker", "run", "-i", img_name],
        capture_output=True,
        text=True,
        input=vars.to_csv(),
    )
    if p.returncode:
        print(f"ERROR predicting resistance status with '{pred_container_tar_path}':")
        print(p.stderr)
        exit(1)
    return bool(p.stdout)

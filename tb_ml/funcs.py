from glob import glob
from uuid import uuid4
from xml.etree.ElementTree import Comment
import pandas as pd
import argparse
from typing import Iterable, Union
import pathogenprofiler as pp
import os
from tqdm import tqdm
import subprocess as sp

def get_cli_args() -> argparse.Namespace:
    pass


def get_positions(AFs: pd.Series) -> list:
    """
    Extracts the genomic positions of interest from the AF series; assumes
    the index entries to have the format: 'POS_REF_ALT'.
    """
    pos = [int(x.split('_')[0]) for x in AFs.index]
    return pos

def get_genotypes(bam_file: str, ref_file: str ,variants: list,subset_bam: bool=True) -> list:
    """
    1. Writes a dummy vcf with the variants 
    2. Subsets the bam file to only include overlapping reads
    3. Uses freebayes to produce calls for those specific postiions. 
    """
    prefix = str(uuid4())
    tmp_vcf_file = f"{prefix}.vcf"
    tmp_bed_file = f"{prefix}.bed"
    tmp_bam_file = f"{prefix}.bam"
    with open(tmp_vcf_file,"w") as O:
        O.write("""##fileformat=VCFv4.2
##reference=/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa
##contig=<ID=Chromosome,length=4411532>
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
""")
        for v in variants:
            pos,ref,alt = v.split("_")
            O.write(f"Chromosome\t{pos}\t.\t{ref}\t{alt}\t.\t.\tAF=1\n")
    if subset_bam:
        sp.call(f"vcf2bed < {tmp_vcf_file} > {tmp_bed_file}",shell=True)
        sp.call(f"samtools view -bML {tmp_bed_file} {bam_file} > {tmp_bam_file}",shell=True)
        bam_file = tmp_bam_file
    rows = []
    for l in tqdm(sp.Popen(
        f"freebayes -f {ref_file} {bam_file} --variant-input {tmp_vcf_file} --only-use-input-alleles  "
        f"| bcftools norm -f {ref_file} -m - "
        "| bcftools query -f '%POS\t%REF\t%ALT\t[%GT\t%DP]\n'",
        shell=True, stdout=sp.PIPE).stdout,total=len(variants)):
        row = l.decode().strip().split()
        varID = "_".join(row[:3])
        rows.append((varID,row[3],row[4]))
    for f in glob(f"{prefix}*"):
        os.remove(f)
    return rows

def run_variant_calling_pipeline(
    bam_file: str, ref_file: str, AFs: pd.Series
) -> pd.DataFrame:
    """
    Runs the variant calling pipeline on a bam file using a
    vector of positions. Fills in with the mean values from AFs
    if a position is missing
    """
    def get_missing_positions(AFs,calls):
        return list(AFs.index.difference(set([c[0] for c in calls])))

    
    calls = []
    while True:
        variants = get_missing_positions(AFs,calls)
        
        calls += get_genotypes(bam_file,ref_file,variants)
        if len(calls)==len(AFs):
            break

    varIDs = []
    genotypes = []
    for var_id,gt,dp in sorted(calls,key=lambda x:list(AFs.index).index(x[0])):
        if gt=="1/1":
            g = 1
        if gt=="0/0":
            g = 0
        else:
            g = AFs.loc[var_id,"mean"]
        varIDs.append(var_id)
        genotypes.append(g)

    genotypes = pd.DataFrame(data={"varID":varIDs,"genotype":genotypes})
    return genotypes


def sanitize_input_dimensions(vars: pd.Series, AFs: pd.Series) -> pd.Series:
    """
    Accepts two Series with indices in the form of `POS_REF_ALT`. Returns a
    new Series safe for passing to the prediction model (i.e. with variants
    not used by the model removed and with AF values for missing variants).
    """
    keep = set(vars.index).intersection(AFs.index)
    missing = set(AFs.index) - set(vars.index)
    new_vars = pd.concat((vars[keep], AFs[missing]))
    # ensure correct order
    new_vars = new_vars[AFs.index]
    return new_vars


def run_prediction_container(container_path: str, vars: pd.Series) -> bool:
    pass


from glob import glob
from uuid import uuid4
from xml.etree.ElementTree import Comment
import pandas as pd
import subprocess
import argparse
from typing import Iterable, Union
import pathogenprofiler as pp
import os
from tqdm import tqdm


def get_cli_args() -> argparse.Namespace:
    pass


def get_positions(AFs: pd.Series) -> list[int]:
    """
    Extracts the genomic positions of interest from the AF series; assumes
    the index entries to have the format: 'POS_REF_ALT'.
    """
    pos = [int(x.split('_')[0]) for x in AFs.index]
    return pos


def get_genotypes(bam_file: str, ref_file: str, variants: list,
                  subset_bam: bool = True) -> pd.Series:
    """
    1. Writes a dummy vcf with the variants
    2. Subsets the bam file to only include overlapping reads
    3. Uses freebayes to produce calls for those specific postiions.
    """
    prefix = str(uuid4())
    tmp_vcf_file = f"{prefix}.vcf"
    tmp_bed_file = f"{prefix}.bed"
    tmp_bam_file = f"{prefix}.bam"
    with open(tmp_vcf_file, "w") as outfile:
        outfile.write('##fileformat=VCFv4.2')
        outfile.write('##reference=/home/jody/refgenome/'
                      'MTB-h37rv_asm19595v2-eg18.fa')
        outfile.write('##contig=<ID=Chromosome,length=4411532>')
        outfile.write('##INFO=<ID=AF,Number=A,Type=Float,Description='
                      '"Estimated allele frequency in the range (0,1]">')
        outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
        for v in variants:
            pos, ref, alt = v.split("_")
            outfile.write(f"Chromosome\t{pos}\t.\t{ref}\t{alt}\t.\t.\tAF=1\n")
    if subset_bam:
        subprocess.run(
            f"vcf2bed < {tmp_vcf_file} > {tmp_bed_file}", shell=True)
        subprocess.run(
            f"samtools view -bML {tmp_bed_file} {bam_file} > {tmp_bam_file}",
            shell=True)
        bam_file = tmp_bam_file
    VC_result = subprocess.run(
        f"freebayes -f {ref_file} {bam_file} --variant-input {tmp_vcf_file} "
        "--only-use-input-alleles  "
        f"| bcftools norm -f {ref_file} -m - "
        "| bcftools query -f '%POS\_%REF\_%ALT\t[%GT\t%DP]\n'",
        shell=True, capture_output=True, text=True).stdout.strip()
    vars = pd.Series(VC_result.split('\n')).str.split('\t', expand=True)
    vars.columns = ['varID', 'GT', 'DP']
    vars = vars.set_index('varID')
    vars['GT'] = vars['GT'].apply(lambda x: int(x.split('/')[0]) if
                                  x != '.' else pd.NA)
    for f in glob(f"{prefix}*"):
        os.remove(f)
    return vars


def run_variant_calling_pipeline(
        bam_file: str, ref_file: str, AFs: pd.Series) -> pd.DataFrame:
    """
    Runs the variant calling pipeline on a bam file using a
    vector of positions. Fills in with the mean values from AFs
    if a position is missing
    """
    def get_missing_positions(AFs, calls):
        return list(AFs.index.difference(set([c[0] for c in calls])))

    calls: list = []
    while True:
        variants = get_missing_positions(AFs, calls)

        calls += get_genotypes(bam_file, ref_file, variants)
        if len(calls) == len(AFs):
            break

    varIDs = []
    genotypes = []
    for var_id, gt, dp in sorted(calls, key=lambda x: list(AFs.index).index(x[0])):
        if gt == "1/1":
            g = 1
        if gt == "0/0":
            g = 0
        else:
            g = AFs.loc[var_id, "mean"]
        varIDs.append(var_id)
        genotypes.append(g)

    genotypes = pd.DataFrame(data={"varID": varIDs, "genotype": genotypes})
    return genotypes


def TEST_run_variant_calling_pipeline(mock_vars_fname: str,
                                      pos: Iterable[int]) -> pd.Series:
    """
    Runs the variant calling pipeline on the raw reads of the sample using a
    vector of positions.
    """
    vars = pd.read_csv(mock_vars_fname, index_col=0).squeeze()
    return vars


def sanitize_input_dimensions(vars: pd.Series, AFs: pd.Series) -> pd.Series:
    """
    Accepts two Series with indices in the form of `POS_REF_ALT`. Returns a
    new Series safe for passing to the prediction model(i.e. with variants
    not used by the model removed and with AF values for missing variants).
    """
    keep = list(set(vars.index).intersection(AFs.index))
    missing = list(set(AFs.index) - set(vars.index))
    new_vars = pd.concat((vars[keep], AFs[missing]))
    # ensure correct order
    new_vars = new_vars[AFs.index]
    return new_vars


def run_prediction_container(pred_container_tar_path: str,
                             vars: pd.Series) -> bool:
    # get the image tag first
    p = subprocess.run(['docker', 'load', '-i', pred_container_tar_path],
                       capture_output=True, text=True)
    img_name = p.stdout.split('Loaded image: ')[1].strip()
    # now run the container got the prediction
    p = subprocess.run(['docker', 'run', '-i', img_name], capture_output=True,
                       text=True, input=vars.to_csv())
    if p.returncode:
        print('ERROR predicting resistance status '
              f'with \'{pred_container_tar_path}\':')
        print(p.stderr)
        exit(1)
    return bool(p.stdout)
from glob import glob
from uuid import uuid4
import pandas as pd
import numpy as np
import subprocess
import argparse
import os


def get_cli_args() -> tuple[str, str, str]:
    parser = argparse.ArgumentParser(
        description="""
        TB-ML: A framework for comparing AMR prediction in M. tuberculosis.
        """
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="Aligned reads for a single sample [required]",
        metavar="FILE",
    )
    parser.add_argument(
        "-p",
        "--prediction-container",
        type=str,
        required=True,
        help="Name of the Docker image to use for prediction [required]",
        metavar="STR",
    )
    parser.add_argument(
        "-vc",
        "--variant-calling-container",
        type=str,
        help="Name of the Docker image to use for variant-calling",
        metavar="STR",
    )

    args = parser.parse_args()
    return args.bam, args.prediction_container, args.variant_calling_container


def dev_test_args():
    import pathlib

    path = f"{pathlib.Path.home()}/git/tb-ml/"
    bam_file = f"{path}/test_data/test.cram"
    af_file = f"{path}/test_data/SM_training_AF.csv"
    vc_container = "vc-test"
    pred_container = "rf-sm-predictor"
    return bam_file, af_file, vc_container, pred_container


def get_positions(
    AFs: pd.Series,
) -> list[int]:
    """
    Extracts the genomic positions of interest from the AF series; assumes
    the index entries to have the format: 'POS_REF_ALT'.
    """
    pos = [int(x.split("_")[0]) for x in AFs.index]
    return pos


def run_VC_container(
    VC_container_img_name: str,
    bam_file: str,
    af_file: str,
    DP_threshold: int = 10,
) -> pd.Series:
    """
    Run a docker container containing a variant-calling pipeline and return the found
    variants in a `pd.Series` with AF values replacing non-calls.
    """
    # run the container (the bind volume needs absolute paths)
    af_path = af_file if os.path.isabs(af_file) else f"{os.getcwd()}/{af_file}"
    bam_path = bam_file if os.path.isabs(bam_file) else f"{os.getcwd()}/{bam_file}"
    p = subprocess.run(
        [
            "docker",
            "run",
            "--mount",
            f"type=bind,source={af_path},target=/data/AFs.csv",
            "--mount",
            f"type=bind,source={bam_path},target=/data/aligned_reads",
            VC_container_img_name,
        ],
        capture_output=True,
        text=True,
    )
    if p.returncode:
        print(f"ERROR running variant-calling pipeline with '{VC_container_img_name}':")
        print(p.stderr)
        exit(1)
    # reformat the result
    variants = pd.Series(p.stdout.strip().split("\n")).str.split("\t", expand=True)
    variants.columns = ["varID", "GT", "DP"]
    variants = variants.set_index("varID")
    variants["GT"] = variants["GT"].apply(lambda x: x[0])
    variants = variants.replace(".", np.nan).astype(float)
    # declare variants with DP < threshold as non-calls
    variants.loc[variants["DP"] < DP_threshold, "GT"] = np.nan
    return variants["GT"]


def process_variants(
    variants: pd.Series,
    AFs: pd.Series,
) -> pd.Series:
    """
    Makes sure `variants` and `AFs` have the same dimensions by dropping entries not
    present in `AFs` and using values from `AFs` for entries missing in `variants`.
    """
    new_vars = variants.copy()
    # replace noncalls with the corresponding AF values
    new_vars.fillna(AFs, inplace=True)
    # add the allele frequencies for variants not found in the variant
    # calling pipeline to ensure matching dimensions for the prediction model
    new_vars = pd.concat((new_vars, AFs[AFs.index.difference(new_vars.index)]))
    # make sure the order is as expected by the model
    new_vars = new_vars[AFs.index]
    return new_vars


def run_prediction_container(
    pred_container_img_name: str,
    variants: pd.Series,
) -> bool:
    """
    Run a docker container containing a model for predicting AMR resistance from a
    `pd.Series` of genotypes.
    """
    # run the container to get the prediction
    p = subprocess.run(
        ["docker", "run", "-i", pred_container_img_name],
        capture_output=True,
        text=True,
        input=variants.to_csv(),
    )
    if p.returncode:
        print(f"ERROR predicting resistance status with '{pred_container_img_name}':")
        print(p.stderr)
        exit(1)
    return bool(p.stdout)

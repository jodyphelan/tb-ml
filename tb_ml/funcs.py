from glob import glob
from uuid import uuid4
import pandas as pd
import numpy as np
from typing import Union
import subprocess
import argparse
import io
import os


def get_cli_args() -> tuple[str, str, str, Union[str, None]]:
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
        "-v",
        "--variant-calling-container",
        type=str,
        required=True,
        help="Name of the Docker image to use for variant-calling [required]",
        metavar="STR",
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
        "--variants-filename",
        type=str,
        help="Write the called variants to this file before prediction [optional]",
        metavar="STR",
    )
    args = parser.parse_args()
    return (
        args.bam,
        args.prediction_container,
        args.variant_calling_container,
        args.variants_filename,
    )


def dev_test_args() -> tuple[str, str, str, str]:
    import pathlib

    path = pathlib.Path(__file__).resolve().parent.parent
    print(path)
    bam_file = f"{path}/test_data/test.cram"
    vc_container = "vc-test"
    pred_container = "rf-sm-predictor"
    write_vars_fname = f"{path}/test_data/test_variants.vcf"
    return bam_file, vc_container, pred_container, write_vars_fname


def run_VC_container(
    VC_container_img_name: str,
    bam_file: str,
    target_vars: pd.Series,
    DP_threshold: int = 10,
) -> pd.Series:
    """
    Run a docker container containing a variant-calling pipeline and return the found
    variants in a `pd.Series`. The container expects the target variants in the format
    'POS,REF,ALT' in STDIN in order to make sure that they are covered in the results.
    """
    # bring the target variants into the right format
    target_vars_str = target_vars.reset_index()[  # type:ignore
        ["POS", "REF", "ALT"]
    ].to_csv(header=False, index=False)
    # run the container (the bind volume needs absolute paths)
    bam_path = bam_file if os.path.isabs(bam_file) else f"{os.getcwd()}/{bam_file}"
    p = subprocess.run(
        [
            "docker",
            "run",
            "-i",
            "--mount",
            f"type=bind,source={bam_path},target=/data/aligned_reads",
            VC_container_img_name,
        ],
        capture_output=True,
        text=True,
        input=target_vars_str,
    )
    if p.returncode:
        print(f"ERROR running variant-calling pipeline with '{VC_container_img_name}':")
        print(p.stderr)
        exit(1)
    # reformat the result
    variants = pd.read_csv(io.StringIO(p.stdout), header=None)
    variants.columns = ["POS", "REF", "ALT", "GT", "DP"]
    variants = variants.set_index(["POS", "REF", "ALT"])
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


def get_target_vars_from_prediction_container(
    pred_container_img_name: str,
) -> pd.Series:
    """
    Run the prediction container with "get_AFs" as argument in order to get the allele
    frequencies of the training dataset. These can then be used for targeted variant
    calling.
    """
    p = subprocess.run(
        ["docker", "run", pred_container_img_name, "get_target_vars"],
        capture_output=True,
        text=True,
    )
    if p.returncode:
        print(f"ERROR extracting allele frequencies from '{pred_container_img_name}':")
        print(p.stderr)
        exit(1)
    target_vars = (
        pd.read_csv(io.StringIO(p.stdout)).set_index(["POS", "REF", "ALT"]).squeeze()
    )
    return target_vars


def run_prediction_container(
    pred_container_img_name: str,
    variants: pd.Series,
) -> float:
    """
    Run a docker container containing a model for predicting AMR resistance from a
    `pd.Series` of genotypes. The result is the probability of resistance.
    """
    # run the container to get the prediction
    p = subprocess.run(
        ["docker", "run", "-i", pred_container_img_name, "predict"],
        capture_output=True,
        text=True,
        input=variants.to_csv(),
    )
    if p.returncode:
        print(f"ERROR predicting resistance status with '{pred_container_img_name}':")
        print(p.stderr)
        exit(1)
    return float(p.stdout.strip())

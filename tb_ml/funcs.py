import pandas as pd
from typing import Optional
import subprocess
import argparse
import io
import os
import sys


DEFAULT_VC_CONTAINER = "julibeg/tb-ml-variant-calling:0.1.1"
DEFAULT_PRED_CONTAINER = "julibeg/tb-ml-streptomycin-rf-predictor:0.1.1"


def get_cli_args() -> tuple[str, str, str, Optional[str], Optional[str]]:
    parser = argparse.ArgumentParser(
        description="""
        TB-ML: A framework for comparing AMR prediction in M. tuberculosis.
        """,
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
        default=DEFAULT_VC_CONTAINER,
        help='Name of the Docker image for variant-calling (default: "%(default)s")',
        metavar="STR",
    )
    parser.add_argument(
        "-p",
        "--prediction-container",
        type=str,
        default=DEFAULT_PRED_CONTAINER,
        help='Name of the Docker image for prediction (default: "%(default)s")',
        metavar="STR",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Write output to this file instead of STDOUT",
        metavar="STR",
    )
    parser.add_argument(
        "--variants-filename",
        type=str,
        help="Write the called variants to this file before prediction",
        metavar="STR",
    )
    args = parser.parse_args()
    return (
        args.bam,
        args.variant_calling_container,
        args.prediction_container,
        args.output,
        args.variants_filename,
    )


def get_prediction(
    bam_file: str,
    vc_container: str,
    pred_container: str,
    write_vars_fname: Optional[str] = None,
) -> pd.Series:
    """
    Main function; uses two Docker containers (one containing a variant-calling pipeline
    and one with a AMR prediction model) to predict AMR resistance from the aligned
    reads of a sample
    (other input types will be added in the future).
    """
    # get the allele frequencies of the training dataset from the prediction container
    target_vars_AF: pd.Series = get_target_vars_from_prediction_container(
        pred_container
    )
    # run the variant calling container and provide the target variants so that it can
    # make sure that they are covered in the output
    vc_stats: pd.Series
    variants: pd.Series
    vc_stats, variants = run_VC_container(bam_file, vc_container, target_vars_AF)
    if write_vars_fname is not None:
        variants.to_csv(write_vars_fname)
    # generate Series for final report
    res = pd.Series(dtype=object)
    res.index.name = "parameter"
    res.name = "value"
    res["file"] = get_absolute_path(bam_file)
    res["vc_container"] = vc_container
    res["pred_container"] = pred_container
    res = pd.concat((res, vc_stats))
    # predict
    res["resistance_probability"] = run_prediction_container(pred_container, variants)
    res["resistance_status"] = "S" if res["resistance_probability"] < 0.5 else "R"
    return res


def get_absolute_path(path: str) -> str:
    """
    Gets the absolute path and is aware of `$HOME` and `~/`
    """
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def run_VC_container(
    bam_file: str,
    VC_container_img_name: str,
    target_vars_AF: pd.Series,
    DP_threshold: int = 10,
) -> tuple[pd.Series, pd.Series]:
    """
    Run a docker container containing a variant-calling pipeline and return the found
    variants in a `pd.Series`. The container expects the target variants in the format
    'POS,REF,ALT' in STDIN in order to make sure that they are covered in the results.
    """
    # run the container (the bind volume needs absolute paths)
    bam_path = get_absolute_path(bam_file)
    p = subprocess.run(
        [
            "docker",
            "run",
            "-i",
            "--mount",
            f"type=bind,source={bam_path},target=/data/aligned_reads",
            VC_container_img_name,
        ],
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        text=True,
        input=target_vars_AF.to_csv(),
    )
    if p.returncode:
        print(f"ERROR running variant-calling pipeline with '{VC_container_img_name}'")
        sys.exit(p.returncode)
    # the first few lines of the result will start with '#' and hold some basic stats
    # about the variant calling process
    stats_lines = []
    for line in p.stdout.split("\n"):
        if line.startswith("#"):
            stats_lines.append(line[1:])
        else:
            break
    stats: pd.Series = pd.read_csv(
        io.StringIO("\n".join(stats_lines)), index_col=0
    ).squeeze()
    variants: pd.Series = pd.read_csv(
        io.StringIO(p.stdout),
        index_col=["POS", "REF", "ALT"],
        comment="#",
    ).squeeze()
    return stats, variants


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
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        text=True,
    )
    if p.returncode:
        print(f"ERROR extracting allele frequencies from '{pred_container_img_name}'")
        sys.exit(p.returncode)
    target_vars: pd.Series = (
        pd.read_csv(io.StringIO(p.stdout)).set_index(["POS", "REF", "ALT"]).squeeze()
    )
    # make sure we actually got a Series
    assert isinstance(target_vars, pd.Series), (
        f"Prediction container {pred_container_img_name} returned "
        "ill-formatted target variants."
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
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        text=True,
        input=variants.to_csv(),
    )
    if p.returncode:
        print(f"ERROR predicting resistance status with '{pred_container_img_name}'")
        sys.exit(p.returncode)
    return float(p.stdout.strip())

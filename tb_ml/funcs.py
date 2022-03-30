from xml.etree.ElementTree import Comment
import pandas as pd
import subprocess
import argparse
from typing import Iterable, Union


def get_cli_args() -> argparse.Namespace:
    pass


def get_positions(AFs: pd.Series) -> list[int]:
    """
    Extracts the genomic positions of interest from the AF series; assumes
    the index entries to have the format: 'POS_REF_ALT'.
    """
    pos = [int(x.split('_')[0]) for x in AFs.index]
    return pos


def run_variant_calling_pipeline(seqfile: str,
pos: Iterable[int]) -> pd.Series:
    """
    Runs the variant calling pipeline on the raw reads of the sample using a
    vector of positions.
    """
    pass


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
    new Series safe for passing to the prediction model (i.e. with variants
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
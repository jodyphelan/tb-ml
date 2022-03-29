from xml.etree.ElementTree import Comment
import pandas as pd
import argparse
from typing import Iterable, Union


def get_cli_args() -> argparse.Namespace:
    pass


def get_positions(AFs: pd.Series) -> list:
    """
    Extracts the genomic positions of interest from the AF series; assumes
    the index entries to have the format: 'POS_REF_ALT'.
    """
    pos = [int(x.split('_')[0]) for x in AFs.index]
    return pos


def run_variant_calling_pipeline(seqs_fname: str, pos: Iterable) -> None:
    """
    Runs the variant calling pipeline on the raw reads of the sample using a
    vector of positions.
    """
    pass


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


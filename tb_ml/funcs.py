import pandas as pd
from typing import Optional
import argpass
import io

from . import util

DEFAULT_VC_CONTAINER = "julibeg/tb-ml-variant-calling"
DEFAULT_PRED_CONTAINER = "julibeg/tb-ml-streptomycin-rf-predictor"


class VariantCallingContainer(util.DockerImage):
    def run_vc_pipeline(
        self,
        bam_file: str,
        target_vars_AF: pd.Series,
        extra_args: Optional[list[str]] = None,
    ) -> tuple[pd.Series, pd.Series]:
        bam_path = util.get_absolute_path(bam_file)
        result = self.run(
            docker_args=[
                "--mount",
                f"type=bind,source={bam_path},target=/data/aligned_reads",
            ],
            extra_args=extra_args,
            input=target_vars_AF.to_csv(),
        )
        # the first few lines of the result will start with '#' and hold some basic
        # stats about the variant calling process
        stats_lines = []
        for line in result.split("\n"):
            if line.startswith("#"):
                stats_lines.append(line[1:])
            else:
                break
        stats: pd.Series = pd.read_csv(
            io.StringIO("\n".join(stats_lines)), index_col=0
        ).squeeze()
        variants: pd.Series = pd.read_csv(
            io.StringIO(result),
            index_col=["POS", "REF", "ALT"],
            comment="#",
        ).squeeze()
        return stats, variants


class PredictionContainer(util.DockerImage):
    def get_target_variants_and_AFs(self) -> pd.Series:
        target_vars: pd.Series = (
            pd.read_csv(io.StringIO(self.run(extra_args=["get_target_vars"])))
            .set_index(["POS", "REF", "ALT"])
            .squeeze()
        )
        return target_vars

    def predict(self, variants: pd.Series) -> float:
        return float(self.run(input=variants.to_csv(), extra_args=["predict"]).strip())


def get_cli_args() -> tuple[
    str,
    str,
    str,
    Optional[str],
    Optional[str],
    Optional[list[str]],
    Optional[list[str]],
]:
    parser = argpass.ArgumentParser(
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
    parser.add_argument(
        "--variant-calling-args",
        type=str,
        nargs=argpass.NargsOption.COLLECT_UNTIL_NEXT_KNOWN,
        help="Extra argument(s) to pass on to the VC container",
        metavar="STR",
    )
    parser.add_argument(
        "--prediction-args",
        type=str,
        nargs=argpass.NargsOption.COLLECT_UNTIL_NEXT_KNOWN,
        help="Extra argument(s) to pass on to the prediction container",
        metavar="STR",
    )
    args = parser.parse_args()
    return (
        args.bam,
        args.variant_calling_container,
        args.prediction_container,
        args.output,
        args.variants_filename,
        args.variant_calling_args,
        args.prediction_args,
    )


def get_prediction(
    bam_file: str,
    vc_img_name: str,
    pred_img_name: str,
    write_vars_fname: Optional[str] = None,
    vc_extra_args: Optional[list[str]] = None,
    pred_extra_args: Optional[list[str]] = None,
) -> pd.Series:
    """
    Main function; uses two Docker containers (one containing a variant-calling pipeline
    and one with a AMR prediction model) to predict AMR resistance from the aligned
    reads of a sample
    (other input types will be added in the future).
    """
    # pull / update the images and initialise the Docker objects
    vc_container: VariantCallingContainer = VariantCallingContainer(vc_img_name)
    pred_container: PredictionContainer = PredictionContainer(pred_img_name)
    # get the allele frequencies of the training dataset from the prediction container
    target_vars_AF: pd.Series = pred_container.get_target_variants_and_AFs()
    # run the variant calling container and provide the target variants so that it can
    # make sure that they are covered in the output
    vc_stats: pd.Series
    variants: pd.Series
    vc_stats, variants = vc_container.run_vc_pipeline(
        bam_file, target_vars_AF, extra_args=vc_extra_args
    )
    if write_vars_fname is not None:
        variants.to_csv(write_vars_fname)
    # generate Series for final report
    res = pd.Series(dtype=object)
    res.index.name = "parameter"
    res.name = "value"
    res["file"] = util.get_absolute_path(bam_file)
    res["vc_container"] = vc_img_name
    res["pred_container"] = pred_img_name
    res = pd.concat((res, vc_stats))
    # predict
    res["resistance_probability"] = pred_container.predict(variants)
    res["resistance_status"] = "S" if res["resistance_probability"] < 0.5 else "R"
    return res

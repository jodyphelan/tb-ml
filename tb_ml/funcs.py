import pandas as pd
from typing import List, Optional
import argpass
import io

from . import util


DEFAULT_VC_CONTAINER = "julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0"
DEFAULT_PRED_CONTAINER = "julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.1.0"


class VariantCallingContainer(util.DockerImage):
    def run_vc_pipeline(
        self,
        bam_file: str,
        target_vars_AF: pd.Series,
        extra_args: Optional[List[str]] = None,
    ) -> tuple[pd.Series, pd.Series]:
        # docker needs absolute paths for mounts
        bam_path = util.get_absolute_path(bam_file)
        # we need to write the target vars to a temporary file
        with util.temp_file() as tmp_path:
            target_vars_AF.to_csv(tmp_path)
            # define mount points in container
            container_bam_path = "/data/aligned_reads"
            container_target_vars_path = "/data/target_vars_AF"
            # define the arguments to be passed to the entrypoint in the container
            extra_args = [] if extra_args is None else extra_args
            extra_args += ["-b", container_bam_path, "-t", container_target_vars_path]
            result = self.run(
                docker_args=[
                    "--mount",
                    f"type=bind,source={bam_path},target={container_bam_path}",
                    "--mount",
                    f"type=bind,source={tmp_path},target={container_target_vars_path}",
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

    def predict(
        self, variants: pd.Series, extra_args: Optional[List[str]] = None
    ) -> float:
        extra_args = [] if extra_args is None else extra_args
        with util.temp_file() as tmp_path:
            # define the arguments to be passed to the entrypoint in the container
            variants.to_csv(tmp_path)
            container_vars_path = "/data/variants.csv"
            extra_args += ["predict", container_vars_path]
            result = self.run(
                extra_args=extra_args,
                docker_args=[
                    "--mount",
                    f"type=bind,source={tmp_path},target={container_vars_path}",
                ],
            ).strip()
        return float(result)


def get_cli_args() -> tuple[str, List[str], str, List[str]]:
    parser = argpass.ArgumentParser(
        description="""
        TB-ML: A framework for comparing AMR prediction in M. tuberculosis.
        """,
    )
    parser.add_argument(
        "--preprocessing-container",
        type=str,
        required=True,
        help='Name of the Docker image for pre-processing [required]',
        metavar="STR",
        dest="preproc_container",
    )
    parser.add_argument(
        "--preprocessing-args",
        type=str,
        required=True,
        help="String with extra argument(s) to pass on to the pre-processing container",
        metavar="STR",
        dest="preproc_args",
    )
    parser.add_argument(
        "--prediction-container",
        type=str,
        required=True,
        help='Name of the Docker image for prediction [required]',
        metavar="STR",
        dest="pred_container",
    )
    parser.add_argument(
        "--prediction-args",
        type=str,
        required=True,
        help="String with extra argument(s) to pass on to the prediction container",
        metavar="STR",
        dest="pred_args",
    )
    args = parser.parse_args()
    preproc_args = [] if args.preproc_args is None else args.preproc_args.split()
    pred_args = [] if args.pred_args is None else args.pred_args.split()
    return (
        args.preproc_container,
        preproc_args,
        args.pred_container,
        pred_args,
    )


def get_prediction(
    bam_file: str,
    vc_img_name: str,
    pred_img_name: str,
    write_vars_fname: Optional[str] = None,
    vc_extra_args: Optional[List[str]] = None,
    pred_extra_args: Optional[List[str]] = None,
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
    res["resistance_probability"] = pred_container.predict(variants, pred_extra_args)
    res["resistance_status"] = "S" if res["resistance_probability"] < 0.5 else "R"
    return res

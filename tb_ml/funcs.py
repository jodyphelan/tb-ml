import pandas as pd
from typing import List, Dict, Tuple, Optional
import argparse
import tempfile
import pathlib
import os
import sys
import time
import io

from . import util


DEFAULT_VC_CONTAINER = "julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0"
DEFAULT_PRED_CONTAINER = "julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.1.0"

output_comment_lines: List[str] = []


def get_cli_args() -> List[Tuple[str, List[str]]]:
    # We will use argparse only to provide the help string, but parse the command line
    # arguments ourselves.
    parser = argparse.ArgumentParser(
        description="""
        TB-ML: A framework for comparing AMR prediction in M. tuberculosis. Provide
        Docker image names and arguments like so: --container CONTR_NAME_1 ARG_1 ARG_2
        --container CONTR_NAME_2 ARG_3 --container CONTR_NAME_3 ARG_4 ARG_5 ...
        """,
    )
    parser.add_argument(
        "--container",
        type=str,
        required=True,
        action="append",
        nargs="+",
        help="Name of Docker image and corresponding extra arguments [required]",
        metavar="STR",
    )
    args, other_args = parser.parse_args()
    # now parse the command line args --> split into lists (with the container name at
    # the beginning followed by the container args) at every "--container" argument
    cont_args_lists: List[List[str]] = []
    args_list: List[str] = []
    for arg in sys.argv[1:]:
        if arg == "--container":
            if args_list:
                cont_args_lists.append(args_list)
                args_list = []
            continue
        args_list.append(arg)
    cont_args_lists.append(args_list)
    containers_and_args = [(lst[0], lst[1:]) for lst in cont_args_lists]
    return containers_and_args


def get_prediction(containers_and_args: List[Tuple[str, List[str]]]) -> pd.Series:
    """
    TODO: add docstring
    """
    # run all container commands in a temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        print(tmp_dir)
        for img_name, container_args in containers_and_args:
            container = util.DockerImage(img_name)
            docker_args = ["-v", f"{tmp_dir}:/data"]
            for i, arg in enumerate(container_args):
                # Check if an argument is a path. If it is, it can either be an input
                # file / directory outside the temporary directory or an intermediate
                # file / directory which has been generated inside the temporary
                # directory. In the first case, we will make it available to the
                # container via a bind mount. In the second case, we don't need to do
                # anything since the temporary directory will be mounted as a volume.
                if os.path.exists(arg):
                    # The path exists, now check if it is inside the temporary
                    # directory. To achieve this, prepend the basename with the path to
                    # the temp dir and check if it exists.
                    abs_path_tmpdir = pathlib.Path(tmp_dir, pathlib.Path(arg).name)
                    if not abs_path_tmpdir.exists():
                        # define the bind mount
                        docker_args += [
                            "--mount",
                            (
                                f"type=bind,source={util.get_absolute_path(arg)},"
                                f"target=/data/{arg}"
                            ),
                        ]
                        # replace the original path with the path to the bind mount
                        container_args[i] = f"/data/{arg}"
            print(img_name, container_args)
            container.run(docker_args=docker_args, extra_args=container_args)
        time.sleep(10000)

    # generate Series for final report
    # res = pd.Series(dtype=object)
    # res.index.name = "parameter"
    # res.name = "value"
    # res["file"] = util.get_absolute_path(bam_file)
    # res["vc_container"] = vc_img_name
    # res["pred_container"] = pred_img_name
    # res = pd.concat((res, vc_stats))
    # # predict
    # res["resistance_probability"] = pred_container.predict(variants, pred_extra_args)
    # res["resistance_status"] = "S" if res["resistance_probability"] < 0.5 else "R"
    # return res

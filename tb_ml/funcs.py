import pandas as pd
from typing import List, Tuple, Optional
import argparse
import tempfile
import pathlib
import os

from . import util


DEFAULT_VC_CONTAINER = "julibeg/tb-ml-freebayes-vc-from-cram:v0.1.0"
DEFAULT_PRED_CONTAINER = "julibeg/tb-ml-simple-rf-predictor-streptomycin:v0.1.0"


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
        nargs=2,
        help="Name of Docker image and corresponding extra arguments [required]",
        metavar="STR",
    )
    args = parser.parse_args()
    container_args_list = [(lst[0], lst[1].strip().split()) for lst in args.container]
    return container_args_list


def get_output_file(args: List[str]) -> Optional[str]:
    """
    Get the argument following a '-o' or '--output' flag in a list of command line
    arguments.
    """
    for i, arg in enumerate(args):
        if arg in ("-o", "--output"):
            return args[i + 1]
    return None


def run_containers(containers_and_args: List[Tuple[str, List[str]]]) -> pd.Series:
    """
    TODO: add docstring
    """
    # run all container commands in a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        print("# tb-ml was run with the following containers + arguments:")
        process_stdout_list: List[str] = []
        for img_name, container_args in containers_and_args:
            print(f"# Container: {img_name}; Args: \"{' '.join(container_args)}\"")
            container = util.DockerImage(img_name)
            docker_args = ["--rm", "-v", f"{temp_dir}:/data"]
            for i, arg in enumerate(container_args):
                # Check if an argument is a path. If it is, it can either be an input
                # file / directory outside the temporary directory or an intermediate
                # file / directory which has been generated inside the temporary
                # directory. In the first case, we will make it available to the
                # container via an extra volume. In the second case, we don't need to do
                # anything since the temporary directory will be mounted as a volume
                # anyway.
                if os.path.exists(arg):
                    # The path exists, now check if it is inside the temporary
                    # directory. To achieve this, prepend the basename with the path to
                    # the temp dir and check if it exists.
                    abs_path_tmpdir = pathlib.Path(temp_dir, pathlib.Path(arg).name)
                    if not abs_path_tmpdir.exists():
                        # add the volume to the docker arguments
                        source_path = util.get_absolute_path(arg)
                        target_path = f"/extra_mounts/{arg}"
                        docker_args += ["-v", f"{source_path}:{target_path}"]
                        # replace the original path with the path to the volume
                        container_args[i] = target_path
            # print(f"Docker args: {' '.join(docker_args)}")
            output = container.run(
                docker_args=docker_args, extra_args=container_args
            )
            process_stdout_list.append(output)
        # concatenate the output of the containers and print to STDOUT to finish the
        # report
        print("# The results are reported below:")
        print("".join(process_stdout_list), end="")

from typing import List, Tuple, Optional
import argparse
import tempfile
import pathlib
import os

from . import util


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


def run_containers(containers_and_args: List[Tuple[str, List[str]]]) -> None:
    """
    Main function running the containers and corresponding arguments in succession. Will
    print a header line for each container. Any output written by a container to STDOUT
    is collected and printed at the end.
    """
    # run all container commands in a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # print first header line
        print("# tb-ml was run with the following containers + arguments:")
        # define list to collect all output (written to STDOUT) from the containers
        process_stdout_list: List[str] = []
        # run all containers
        for img_name, container_args in containers_and_args:
            # add the name of the container and the arguments to the header
            print(f"# Container: {img_name}; Args: \"{' '.join(container_args)}\"")
            # initialise the Docker wrapper
            container = util.DockerImage(img_name)
            # we will mount the temp dir to the container at `/data`
            docker_args = ["--rm", "-v", f"{temp_dir}:/data"]
            # For each argument, check if it is a path. If it is, it can either be an
            # original input file from outside the temporary directory or an
            # intermediate file which has been generated inside the temporary directory.
            # In the first case, we will make it available to the container via an extra
            # volume. In the second case, we don't need to do anything since the
            # temporary directory will be mounted at `/data` anyway.
            for i, arg in enumerate(container_args):
                if os.path.exists(util.resolve_path(arg)):
                    # The path exists, now check if it is inside the temporary
                    # directory. To achieve this, prepend the basename with the path of
                    # the temp dir and check if this exists.
                    abs_path_tmpdir = pathlib.Path(temp_dir, pathlib.Path(arg).name)
                    if not abs_path_tmpdir.exists():
                        # the argument is a file that's outside of the temp dir --> add
                        # extra volume to the docker arguments
                        source_path = util.get_absolute_path(arg)
                        target_path = f"/extra_mounts/{arg}"
                        docker_args += ["-v", f"{source_path}:{target_path}"]
                        # replace the original path in the argument list with the path
                        # to the volume
                        container_args[i] = target_path
            # now run the container and add any output to the list
            output = container.run(docker_args=docker_args, extra_args=container_args)
            process_stdout_list.append(output)
        # concatenate the output of the containers and print to finish the report
        print("# The results are reported below:")
        print("".join(process_stdout_list), end="")

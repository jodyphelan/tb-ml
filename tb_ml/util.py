import subprocess
from typing import Optional
import os


def get_absolute_path(path: str) -> str:
    """
    Gets the absolute path and is aware of `$HOME` and `~/`
    """
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


class DockerError(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__(msg)


class DockerImage:
    def __init__(self, img_name: str) -> None:
        # check if the image name has a tag. Add 'latest' if not
        self.img_name: str = (f"{img_name}:latest") if ":" not in img_name else img_name

    def exec_cmd(
        self,
        docker_args: list[str],
        extra_args: Optional[list[str]] = None,
        input: Optional[str] = None,
        error_msg: Optional[str] = None,
    ) -> str:
        # replace the default `None` with an empty list
        extra_args = [] if extra_args is None else extra_args
        # the container is supposed to read from STDIN --> make sure `-i` is in the
        # arguments
        if input is not None and "-i" not in docker_args:
            docker_args = docker_args + ["-i"]
        p = subprocess.run(
            [
                "docker",
                *docker_args,
                self.img_name,
                *extra_args,
            ],
            input=input,
            capture_output=True,
            text=True,
        )
        if p.returncode:
            raise DockerError(f"{error_msg}:\n{p.stderr}")
        return p.stdout

    def run(
        self,
        docker_args: Optional[list[str]] = None,
        extra_args: Optional[list[str]] = None,
        input: Optional[str] = None,
        error_msg: Optional[str] = None,
    ) -> str:
        """
        Thin wrapper around `.exec_cmd()`; executes a `docker run ...` command.
        """
        # make sure the docker command starts with `run`
        docker_args = (
            ["run"]
            if (docker_args is None or len(docker_args) == 0)
            else ["run"] + docker_args
            if docker_args[0] != "run"
            else docker_args
        )
        return self.exec_cmd(docker_args, extra_args, input, error_msg)

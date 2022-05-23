import subprocess
import sys
from typing import Optional
from typing import Iterator
from contextlib import contextmanager
import tempfile
import os


@contextmanager
def temp_file() -> Iterator[str]:
    """
    Context manager for creating and deleting a temporary file (works cross-platform
    as opposed to `tempfile.NamedTemporaryFile`)
    """
    tmp_fd, tmp_path = tempfile.mkstemp()
    os.close(tmp_fd)
    try:
        yield tmp_path
    finally:
        os.remove(tmp_path)


def get_absolute_path(path: str) -> str:
    """
    Gets the absolute path and is aware of `$HOME` and `~/`
    """
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


class DockerError(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__(msg)


class DockerImage:
    def __init__(self, img_name: str, pull: bool = False) -> None:
        # check if the image name has a tag. Add 'latest' if not
        img_name = (f"{img_name}:latest") if ":" not in img_name else img_name
        self.img_name: str = img_name
        # pull from the repo, if the image is already present on the host and
        # up to date, `docker pull` won't do anything.
        if pull:
            self.pull()

    def pull(self) -> None:
        """
        Pull the image from the repo. If the local version is already up to date,
        `docker pull` won't do anything.
        """
        output: str = self.exec_cmd(
            ["pull"],
            error_msg=f'Failed pulling image "{self.img_name}"',
        )
        if "Status: Image is up to date" in output:
            print(f'Local image "{self.img_name}" is up to date', file=sys.stderr)

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
        Thin wrapper around `.exec_cmd()` running a command with `docker run ...`
        """
        # make sure the docker command starts with `run`
        docker_args = ["run"] if docker_args is None else ["run"] + docker_args
        return self.exec_cmd(docker_args, extra_args, input, error_msg)

import tb_ml


def main() -> None:
    # process the CLI args to get Docker image names and arguments
    containers_and_args = tb_ml.get_cli_args()
    # run the analysis
    tb_ml.run_containers(containers_and_args)


if __name__ == "__main__":
    main()

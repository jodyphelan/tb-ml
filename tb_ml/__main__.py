import tb_ml


def main() -> None:
    # process CLI args
    # (preproc_container, preproc_args, pred_container, pred_args) = tb_ml.get_cli_args()
    containers_and_args = tb_ml.get_cli_args()
    # run analysis
    res = tb_ml.get_prediction(containers_and_args)
    # print(res.to_csv(), end="")


if __name__ == "__main__":
    main()

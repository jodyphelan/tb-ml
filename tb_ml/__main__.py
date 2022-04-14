import tb_ml


def main() -> None:
    print(tb_ml.get_cli_args())
    # process CLI args
    # (
    #     bam_file,
    #     vc_container,
    #     pred_container,
    #     outfile,
    #     write_vars_fname,
    # ) = tb_ml.get_cli_args()
    # # run analysis
    # res = tb_ml.get_prediction(bam_file, vc_container, pred_container, write_vars_fname)
    # if outfile is not None:
    #     res.to_csv(outfile)
    # else:
    #     print(res.to_csv(), end="")


if __name__ == "__main__":
    main()

import tb_ml


def main() -> None:
    # process CLI args
    (
        bam_file,
        vc_container,
        pred_container,
        outfile,
        write_vars_fname,
        variant_calling_args,
        prediction_args,
    ) = tb_ml.get_cli_args()
    # run analysis
    res = tb_ml.get_prediction(
        bam_file,
        vc_container,
        pred_container,
        write_vars_fname,
        variant_calling_args,
        prediction_args,
    )
    if outfile is not None:
        res.to_csv(outfile)
    else:
        print(res.to_csv(), end="")


if __name__ == "__main__":
    main()

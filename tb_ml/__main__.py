import pandas as pd
import tb_ml


def main() -> None:
    """
    Test what we got so far using hardcoded filenames.
    """
    # process CLI args
    # bam_fname, pred_cont, vc_cont = tb_ml.get_cli_args()
    bam_file, af_file, vc_container, pred_container = tb_ml.dev_test_args()

    AFs: pd.Series = pd.read_csv(af_file, index_col=0).squeeze()
    # run VC pipeline
    variants = tb_ml.run_VC_container(vc_container, bam_file, af_file)
    # process variants to ensure proper dimensions for the prediction model
    variants = tb_ml.process_variants(variants, AFs)
    # predict
    resistance_status: bool = tb_ml.run_prediction_container(pred_container, variants)
    if resistance_status:
        print("resistant")
    else:
        print("not resistant")


if __name__ == "__main__":
    main()

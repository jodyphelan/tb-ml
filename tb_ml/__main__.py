import pandas as pd
import tb_ml


def main() -> None:
    """
    Test what we got so far using hardcoded filenames.
    """
    # process CLI args
    bam_file, pred_container, vc_container = tb_ml.get_cli_args()
    # bam_file, vc_container, pred_container = tb_ml.dev_test_args()

    # get the allele frequencies of the training dataset from the prediction container
    target_vars_AF: pd.Series = tb_ml.get_target_vars_from_prediction_container(
        pred_container
    )
    # run the variant calling container and provide the target variants so that it can 
    # make sure that they are covered in the output
    variants: pd.Series = tb_ml.run_VC_container(vc_container, bam_file, target_vars_AF)
    # process variants to ensure proper dimensions for the prediction model
    variants = tb_ml.process_variants(variants, target_vars_AF)
    # predict
    resistance_status: bool = tb_ml.run_prediction_container(pred_container, variants)
    if resistance_status:
        print("resistant")
    else:
        print("not resistant")


if __name__ == "__main__":
    main()

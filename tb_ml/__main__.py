import pandas as pd
import tb_ml


def main() -> None:
    """
    Test what we got so far using hardcoded filenames.
    """
    # process CLI args
    # bam_file, pred_container, vc_container, write_vars_fname = tb_ml.get_cli_args()
    bam_file, vc_container, pred_container, write_vars_fname = tb_ml.dev_test_args()

    # get the allele frequencies of the training dataset from the prediction container
    target_vars_AF: pd.Series = tb_ml.get_target_vars_from_prediction_container(
        pred_container
    )
    # run the variant calling container and provide the target variants so that it can
    # make sure that they are covered in the output
    variants: pd.Series = tb_ml.run_VC_container(vc_container, bam_file, target_vars_AF)
    # process variants to ensure proper dimensions for the prediction model
    variants_proc = tb_ml.process_variants(variants, target_vars_AF)
    if write_vars_fname is not None:
        variants_proc.to_csv(write_vars_fname)
    res = pd.Series(dtype=float)
    res.name = (
        f"file:{bam_file};vc_container:{vc_container};pred-container:{pred_container}"
    )
    # predict
    res["resistance_probability"] = tb_ml.run_prediction_container(
        pred_container, variants_proc
    )
    # get a few other stats to help interpret the result
    # the number of variants of interest present in the initial genotype array and
    res["shared_variants"] = len(target_vars_AF.index.intersection(variants.index))
    res["dropped_variants"] = len(variants.index.difference(target_vars_AF.index))
    res["variants_set_to_AF"] = len(target_vars_AF.index.difference(variants.index))
    print(res.to_csv(), end="")


if __name__ == "__main__":
    main()
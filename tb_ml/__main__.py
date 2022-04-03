import pandas as pd
import tb_ml


def main() -> None:
    """
    Test what we got so far using hardcoded filenames.
    """
    af_fname = "/home/julsy/git/tb-ml/test_data/SM_training_AF.csv"
    bam_fname = "/home/julsy/git/tb-ml/test_data/test.cram"

    AFs: pd.Series = pd.read_csv(af_fname, index_col=0).squeeze()
    variants = tb_ml.run_VC_container(
        "/home/julsy/git/tb-ml/docker/freebayes_VC_pipeline/VC_container.tar.gz",
        bam_fname,
        af_fname,
    )
    variants = tb_ml.sanitise_variants(variants, AFs)
    pred_container_tar_fname = (
        "/home/julsy/git/tb-ml/docker/simple_RF_predictor/"
        "rf-sm-predictor.docker.tar.gz"
    )
    resistance_status: bool = tb_ml.run_prediction_container(
        pred_container_tar_fname, variants
    )
    if resistance_status:
        print("resistant")
    else:
        print("not resistant")


if __name__ == "__main__":
    main()

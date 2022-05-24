import tb_ml
import pathlib
import pandas as pd
import pandas.testing as pdt


test_data = f"{pathlib.Path(__file__).resolve().parent}/test_data"

VC_CONTAINER_IMG_NAME = tb_ml.DEFAULT_VC_CONTAINER
PRED_CONTAINER_IMG_NAME = tb_ml.DEFAULT_PRED_CONTAINER

bam_file = f"{test_data}/test.cram"


def _get_expected_variants():
    return pd.read_csv(
        f"{test_data}/test_variants.csv", index_col=["POS", "REF", "ALT"]
    ).squeeze()


def test_full():
    """
    Test the full pipeline.
    """
    # get the expected result
    exp_result = pd.read_csv(f"{test_data}/test_result.csv", index_col=0).squeeze()
    # get the actual result
    res = tb_ml.get_prediction(bam_file, VC_CONTAINER_IMG_NAME, PRED_CONTAINER_IMG_NAME)
    # some info in the final report might have changed, but we are only interested in
    # the prediction and VC stats --> drop these rows
    res.drop(["file", "vc_container", "pred_container"], inplace=True)
    # `exp_result` was read into a `Series` with all strings --> change the datatypes
    # to match up with res
    for idx in exp_result.index:
        exp_result[idx] = type(res[idx])(exp_result[idx])
    # now check the results
    pdt.assert_series_equal(exp_result, res)


def test_variants():
    """
    Test the variant calling.
    """
    # get the expected variants
    exp_vars = _get_expected_variants()
    exp_vc_stats = pd.read_csv(f"{test_data}/test_vc_stats.csv", index_col=0).squeeze()
    # initialize the image objects
    vc_container = tb_ml.VariantCallingContainer(VC_CONTAINER_IMG_NAME)
    pred_container = tb_ml.PredictionContainer(PRED_CONTAINER_IMG_NAME)
    # get the target variants and allele frequencies
    target_vars_AF: pd.Series = pred_container.get_target_variants_and_AFs()
    # predict and compare
    vc_stats, variants = vc_container.run_vc_pipeline(bam_file, target_vars_AF)
    pdt.assert_series_equal(exp_vc_stats, vc_stats)
    pdt.assert_series_equal(exp_vars, variants)


def test_prediction():
    """
    Test the prediction
    """
    pred_container = tb_ml.PredictionContainer(PRED_CONTAINER_IMG_NAME)
    variants = _get_expected_variants()
    pred = pred_container.predict(variants)
    assert pred == 0.33

import tb_ml
import pathlib
import pandas as pd
import pandas.testing as pdt


root = pathlib.Path(__file__).resolve().parent.parent.parent

VC_CONTAINER = tb_ml.DEFAULT_VC_CONTAINER
PRED_CONTAINER = tb_ml.DEFAULT_PRED_CONTAINER

bam_file = f"{root}/test_data/test.cram"


def test_full() -> None:
    # get the expected result
    exp_result = pd.read_csv(f"{root}/test_data/test_result.csv", index_col=0).squeeze()
    # get the actual result
    res = tb_ml.get_prediction(bam_file, VC_CONTAINER, PRED_CONTAINER)
    # some info in the final report might have changed, but we are only interested in
    # the prediction and VC stats --> drop these rows
    res.drop(["file", "vc_container", "pred_container"], inplace=True)
    # `exp_result` was read into a `Series` with all strings --> change the datatypes
    # to match up with res
    for idx in exp_result.index:
        exp_result[idx] = type(res[idx])(exp_result[idx])
    # now check the results
    pdt.assert_series_equal(exp_result, res)


def test_variants() -> None:
    # get the expected variants
    exp_vars = pd.read_csv(
        f"{root}/test_data/test_variants.csv", index_col=["POS", "REF", "ALT"]
    ).squeeze()
    exp_vc_stats = pd.read_csv(
        f"{root}/test_data/test_vc_stats.csv", index_col=0
    ).squeeze()
    target_vars_AF = tb_ml.get_target_vars_from_prediction_container(PRED_CONTAINER)
    vc_stats, variants = tb_ml.run_VC_container(bam_file, VC_CONTAINER, target_vars_AF)
    pdt.assert_series_equal(exp_vc_stats, vc_stats)
    pdt.assert_series_equal(exp_vars, variants)

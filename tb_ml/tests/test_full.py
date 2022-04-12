import tb_ml
import pathlib
import pandas as pd
import pandas.testing as pdt


root = pathlib.Path(__file__).resolve().parent.parent.parent

bam_file = f"{root}/test_data/test.cram"


def test_full() -> None:
    # get the expected result
    exp_result = pd.read_csv(f"{root}/test_data/test_result.csv", index_col=0).squeeze()
    # get the actual result
    res = tb_ml.get_prediction(
        bam_file, tb_ml.DEFAULT_VC_CONTAINER, tb_ml.DEFAULT_PRED_CONTAINER
    )
    # `exp_result` was read into a `Series` with
    for idx in exp_result.index:
        exp_result[idx] = type(res[idx])(exp_result[idx])
    # the bam filename in the result `Series` will contain the absolute path, which will
    # be different from the expected results --> drop that row from both `Series`
    pdt.assert_series_equal(exp_result.drop("file"), res.drop("file"))


def test_VC_pipeline() -> None:
    # get the expected variants
    exp_vars = pd.read_csv(
        f"{root}/test_data/test_variants.csv", index_col=["POS", "REF", "ALT"]
    ).squeeze()
    target_vars_AF: pd.Series = tb_ml.get_target_vars_from_prediction_container(
        tb_ml.DEFAULT_PRED_CONTAINER
    )
    variants: pd.Series = tb_ml.run_VC_container(
        bam_file, tb_ml.DEFAULT_VC_CONTAINER, target_vars_AF
    )
    variants_proc: pd.Series = tb_ml.process_variants(variants, target_vars_AF)
    pdt.assert_series_equal(exp_vars, variants_proc)

import tb_ml
import pathlib
import pandas as pd
import pandas.testing as pdt


root = pathlib.Path(__file__).resolve().parent.parent.parent

bam_file = f"{root}/test_data/test.cram"
vc_container = "vc-test"
pred_container = "rf-sm-predictor"


def test_full() -> None:
    # get the expected result
    exp_result = pd.read_csv(f"{root}/test_data/test_result.csv", index_col=0).squeeze()
    # get the actual result
    res = tb_ml.get_prediction(bam_file, vc_container, pred_container)
    # rename the index (the assertion would fail otherwise)
    res.name = "test"
    pdt.assert_series_equal(exp_result, res)


def test_VC_pipeline() -> None:
    # get the expected variants
    exp_vars = pd.read_csv(
        f"{root}/test_data/test_variants.csv", index_col=["POS", "REF", "ALT"]
    ).squeeze()
    target_vars_AF: pd.Series = tb_ml.get_target_vars_from_prediction_container(
        pred_container
    )
    variants: pd.Series = tb_ml.run_VC_container(bam_file, vc_container, target_vars_AF)
    variants_proc: pd.Series = tb_ml.process_variants(variants, target_vars_AF)
    pdt.assert_series_equal(exp_vars, variants_proc)


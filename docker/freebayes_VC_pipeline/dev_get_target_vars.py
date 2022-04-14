import tb_ml


print(
    tb_ml.get_target_vars_from_prediction_container(
        tb_ml.DEFAULT_PRED_CONTAINER
    ).to_csv(),
    end="",
)

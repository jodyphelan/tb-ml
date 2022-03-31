import pandas as pd
from tb_ml import funcs


def main() -> None:
    """
    Test what we got so far using a mock function for step two (the variant
    calling).
    """
    af_fname: str = 'test_data/SM_training_means.csv'
    AFs: pd.Series = pd.read_csv(af_fname, index_col=0).squeeze()
    pos: list[int] = funcs.get_positions(AFs)
    vars = funcs.TEST_run_variant_calling_pipeline(
        'test_data/INH_ERR1035006.csv', pos)
    vars_rdy = funcs.sanitize_input_dimensions(vars, AFs)
    docker_img_fname = \
        'docker/simple_RF_predictor/rf-sm-predictor.docker.tar.gz'
    res_status: bool = funcs.run_prediction_container(docker_img_fname,
                                                      vars_rdy)
    if res_status:
        print('resistant')
    else:
        print('not resistant')


# if __name__ == 'main':
main()
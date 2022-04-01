import pandas as pd
from tb_ml import funcs


def main() -> None:
    """
    Test what we got so far using hardcoded filenames.
    """
    af_fname = 'test_data/SM_training_means.csv'
    bam_fname = 'test_data/test.cram'
    ref_fname = 'test_data/MTB-h37rv_asm19595v2-eg18.fa'

    AFs: pd.Series = pd.read_csv(af_fname, index_col=0).squeeze()
    pos: list[int] = funcs.get_positions(AFs)
    #     'test_data/INH_ERR1035006.csv', pos)
    # vars = get_genotypes()
    variants = funcs.get_genotypes(bam_fname, ref_fname, AFs)
    docker_img_fname = \
        'docker/simple_RF_predictor/rf-sm-predictor.docker.tar.gz'
    resistance_status: bool = funcs.run_prediction_container(docker_img_fname,
                                                             variants)
    if resistance_status:
        print('resistant')
    else:
        print('not resistant')


if __name__ == '__main__':
    main()
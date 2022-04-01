import pandas as pd
import tb_ml


def main() -> None:
    """
    Test what we got so far using hardcoded filenames.
    """
    af_fname = '/home/julian/git/tb-ml/test_data/SM_training_means.csv'
    bam_fname = '/home/julian/git/tb-ml/test_data/test.cram'
    ref_fname = '/home/julian/git/tb-ml/test_data/MTB-h37rv_asm19595v2-eg18.fa'

    AFs: pd.Series = pd.read_csv(af_fname, index_col=0).squeeze()
    variants = tb_ml.get_genotypes(bam_fname, ref_fname, AFs)
    docker_img_fname = ('/home/julian/git/tb-ml/docker/simple_RF_predictor/'
                        'rf-sm-predictor.docker.tar.gz')
    resistance_status: bool = tb_ml.run_prediction_container(
        docker_img_fname, variants)
    if resistance_status:
        print('resistant')
    else:
        print('not resistant')


if __name__ == '__main__':
    main()

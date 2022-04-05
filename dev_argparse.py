import argparse


def get_cli_args() -> tuple[str, str, str]:
    parser = argparse.ArgumentParser(
        description="""
        TB-ML: A framework for comparing AMR prediction in M. tuberculosis.
        """
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="Aligned reads for a single sample [required]",
        metavar="FILE",
    )
    parser.add_argument(
        '-p',
        "--prediction-container",
        type=str,
        required=True,
        help="Name of the Docker image to use for prediction [required]",
        metavar="STR",
        dest='pred_cont',
    )
    parser.add_argument(
        '-vc',
        "--variant-calling-container",
        type=str,
        help="Name of the Docker image to use for variant-calling",
        metavar="STR",
        dest='vc_cont',
    )

    args = parser.parse_args()
    return args.bam, args.pred_cont, args.vc_cont


bam, pred_cont, vc_cont = get_cli_args()
print(bam)
print(pred_cont)
print(vc_cont)
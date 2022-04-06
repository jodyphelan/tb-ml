import pandas as pd
import tb_ml
from typing import Union


def main() -> None:
    # process CLI args
    bam_file, pred_container, vc_container, write_vars_fname = tb_ml.get_cli_args()
    # run analysis
    res = tb_ml.get_prediction(bam_file, vc_container, pred_container, write_vars_fname)
    print(res.to_csv(), end="")


if __name__ == "__main__":
    main()

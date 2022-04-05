import numpy as np
import pandas as pd
import joblib
import sys

"""
Entrypoint for a simple example prediction container with a random forest model fitted
on streptomycin-resistance data. If an argument was passed to it and this argument is
`get_target_vars`, then print the target variants (including AFs in the training set)
to STDOUT so that they can be used to make sure that the variants used by the model are
covered in the variant calling step.
If no argument is passed or the argument is `predict`, read a CSV with genotypes from
STDIN and print the predicted resistance status to STDOUT.
"""

if len(sys.argv) == 2 and sys.argv[1] == "get_target_vars":
    with open("target_vars.csv", "r") as f:
        for line in f:
            sys.stdout.write(line)
elif (len(sys.argv) == 2 and sys.argv[1] == "predict") or len(sys.argv) == 1:
    # load the model
    m = joblib.load("model.pkl")
    # get the prediction input from STDIN
    X = pd.read_csv(sys.stdin, index_col=[0, 1, 2])  # type: ignore
    if X.shape[1] == 1:
        # the DataFrame has only 1 column --> there is only one sample --> it should be
        # the only row --> transpose
        X = X.T
    # predict and print the result
    ypred = m.predict(X)
    print(np.squeeze(ypred))
else:
    sys.exit(
        'ERROR: Provide a single argument ("get_target_vars" or "predict") '
        "or no argument at all."
    )

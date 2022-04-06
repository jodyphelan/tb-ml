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

# get the variants the model has been fitted on
target_vars = pd.read_csv('target_vars.csv', index_col=['POS', 'REF', 'ALT']).squeeze()

if len(sys.argv) == 2 and sys.argv[1] == "get_target_vars":
    sys.stdout.write(target_vars.to_csv())
elif (len(sys.argv) == 2 and sys.argv[1] == "predict") or len(sys.argv) == 1:
    # load the model
    m = joblib.load("model.pkl")
    # get the prediction input from STDIN
    X = pd.read_csv(sys.stdin, index_col=['POS', 'REF', 'ALT'])
    if X.shape[1] == 1:
        # the DataFrame has only 1 column --> there is only one sample --> it should be
        # the only row --> transpose
        X = X.T
    # make sure the variant order is as expected
    X = X.loc[:, target_vars.index]
    # predict the probability for resistance and print the result
    ypred = m.predict_proba(X)[:, 1]
    print(np.squeeze(ypred))
else:
    sys.exit(
        'ERROR: Provide a single argument ("get_target_vars" or "predict") '
        "or no argument at all."
    )

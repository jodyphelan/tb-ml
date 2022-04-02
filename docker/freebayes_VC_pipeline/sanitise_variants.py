import pandas as pd
import sys


af_fname, DP_threshold = sys.argv[1], int(sys.argv[2])

variants = pd.read_csv('new_vars.tsv', sep='\t', index_col=0,
                       header=None)
AFs = pd.read_csv('/data/SM_training_AF.csv', index_col=0, squeeze=True)
variants.columns = ['GT', 'DP']
variants.index.name = 'varID'

# declare variants with DP < threshold as non-calls and replace all
# noncalls with the corresponding AF values
noncalls_idx = [i for i, row in variants.iterrows() if row['GT'] == '.' or
                row['DP'] == '.' or int(row['DP']) < DP_threshold]
variants.loc[noncalls_idx, 'GT'] = AFs[noncalls_idx]
# add the allele frequencies for variants not found in the variant
# calling pipeline to ensure matching dimensions for the prediction model
variants = pd.concat((variants['GT'], AFs[[x for x in AFs.index if x
                                           not in variants.index]]))
# make sure the order is as expected by the model
variants = variants[AFs.index]

# write out to CSV
variants.to_csv('vars_for_prediction.csv')

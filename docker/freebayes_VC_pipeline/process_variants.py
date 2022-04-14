import pandas as pd
import sys

"""
Takes "raw" variants from the variant calling pipeline and makes sure that all variants
used by the predictor are present in the output and that there are no missing data by
replacing non-calls and missing variants with the allele frequencies from the training
dataset (which are passed to the container in STDIN).
Some simple stats (e.g. the number of missing variants etc.) are printed as comments to
STDOUT before the processed variants are printed as CSV in the format `POS,REF,ALT,GT`.
"""

DP_threshold = 10

# read the data
AFs = pd.read_csv("target_vars_AF.csv", index_col=["POS", "REF", "ALT"]).squeeze()
variants = pd.read_csv(sys.stdin, index_col=["POS", "REF", "ALT"])

# get the first char from the `x/x` genotype field and replace non-calls with NA
variants["GT"] = variants["GT"].apply(lambda x: x[0]).replace(".", pd.NA)
# declare variants with DP < threshold also as non-calls
variants.loc[variants["DP"] == ".", "DP"] = -1
variants.loc[variants["DP"].astype(int) < DP_threshold, "GT"] = pd.NA
# collect basic stats on how many variants had to be replaed etc.
stats = pd.Series(dtype=object)
stats["shared_variants"] = len(AFs.index.intersection(variants.index))
stats["dropped_variants"] = len(variants.index.difference(AFs.index))
stats["missing_variants"] = len(AFs.index.difference(variants.index))
stats["noncalls"] = variants["GT"].isna().sum()
stats["variants_set_to_AF"] = stats[["noncalls", "missing_variants"]].sum()
# replace non-calls and missing variants with the corresponding AF values
variants = variants["GT"].fillna(AFs).astype(float)
variants = pd.concat((variants, AFs[AFs.index.difference(variants.index)]))
# make sure the order is as expected by the model
variants = variants[AFs.index]

# write to STDOUT --> first the stats as comments and then the variants
stats_str = stats.to_csv(header=["value"], index_label="parameter")
for line in stats_str.strip().split('\n'):
    print(f'#{line.strip()}')
variants.name = "GT"
print(variants.to_csv(), end="")

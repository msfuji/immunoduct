# Convert immunoduct.gct into an input file of Clustergrammer.
# Masashi Fujita  Aug. 24, 2018

import pandas as pd
from sklearn.preprocessing import MinMaxScaler

def series_to_tuple(s):
    return str(tuple(map(lambda x:"%s: %s" % x, s.items())))

imm_file = snakemake.input["imm"]
ann_file = snakemake.input["ann"]
outfile = snakemake.output["file"]

#
# load GCT file of immune signatures
#
imm = pd.read_csv(imm_file, skiprows=2, sep="\t")

#
# load annotation of samples
#
ann = pd.read_csv(ann_file, sep=None, engine='python')

#
# rename index
#
index = imm[['Name', 'Description']].apply(series_to_tuple, axis=1)
imm.index = index
imm.drop(['Name', 'Description'], axis=1, inplace=True)

#
# match the row order of annotation to the column order of imm
#
ann.index = ann['Sample']
ann = ann.reindex(list(imm.columns))

#
# rename columns
#
columns=ann.apply(series_to_tuple, axis=1)
imm.columns=columns

#
# scaling
#
imm = preprocessing.maxabs_scale(imm, axis=1)

#
# save
#
imm.to_csv(outfile, sep="\t", index_label='')

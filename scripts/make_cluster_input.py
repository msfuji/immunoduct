# Convert immunoduct.gct into an input file of Clustergrammer.
# Masashi Fujita  Aug. 24, 2018

import pandas as pd
import numpy as np
#from sklearn import preprocessing

def series_to_tuple(s):
    return str(tuple(map(lambda x:"%s: %s" % x, s.items())))

def scale(s):
    min = np.min(s)
    max = np.max(s)
    if min > 0 and max > 10*min:
        s = np.log10(s)
        min = np.min(s)
        max = np.max(s)
    elif min == 0 and max > 100:
        s = np.log10(s + 0.01)
        min = np.min(s)
        max = np.max(s)
    s = (s - min)/ (max - min)
    s = 2 * s -1
    return s

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
# make new index
#
index = imm[['Name', 'Description']].apply(series_to_tuple, axis=1)
imm.drop(['Name', 'Description'], axis=1, inplace=True)

#
# match the row order of annotation to the column order of imm
#
ann.index = ann['Sample']
ann = ann.reindex(list(imm.columns))

#
# make new columns
#
columns=ann.apply(series_to_tuple, axis=1)

#
# scaling
#
#imm = preprocessing.minmax_scale(imm, feature_range=(-1, 1), axis=1)
imm = imm.apply(scale, axis=1)
imm.index = index
imm.columns = columns

#imm = pd.DataFrame(imm, index=index, columns=columns)


#
# save
#
imm.to_csv(outfile, sep="\t", index_label='')

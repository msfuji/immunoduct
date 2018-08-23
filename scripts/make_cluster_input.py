import pandas as pd

def read_gct(filename):
    pd.read_csv()
#print(snakemake.input["imm"])
imm_file = "immunoduct.gct"
#print(snakemake.input["ann"])
ann_file = "example/annotation.txt"

#print(snakemake.output["file"])
outfile = "mogemogemoge"

#
# load GCT file of immune signatures
#
imm = pd.read_csv(imm_file, skiprows=2, sep="\t")
print(imm)

#
# load annotation of samples
#
ann = pd.read_csv(ann_file, sep=None)
print(ann)

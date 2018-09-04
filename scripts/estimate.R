source("scripts/common.R")
library(estimate)


args <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]


tmpfile1 <- tempfile()
tmpfile2 <- tempfile()
tmpfile3 <- tempfile()

df <- read_gct(infile)

#
# remove redundant genes
#
df <- df %>% select(-Description) %>% group_by(Name) %>% slice(1) %>% ungroup

#
# save as a matrix
#
cols <- colnames(df)[-1]
cols <- paste(cols, collapse="\t")
cat(cols, "\n", file=tmpfile1)

df %>% write_tsv(tmpfile1, append=T, col_names=F)

#
# filter genes
#
filterCommonGenes(tmpfile1, tmpfile2)

#
# run ESTIMATE
#
estimateScore(tmpfile2, tmpfile3, platform="illumina")

#
# reformat result
#
df <- fread(tmpfile3, skip=2)
df <- df %>% rename(Name=NAME)
df %>% write_gct(outfile)

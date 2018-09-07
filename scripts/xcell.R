source("scripts/common.R")
library(xCell)


args <- commandArgs(trailingOnly=T)
threads <- as.integer(args[1])
infile <- args[2]
outfile <- args[3]

df <- read_gct(infile)

#
# remove redundant genes
#
df <- df %>% group_by(Name) %>% slice(1) %>% ungroup

#
# convert into a matrix
#
mat <- df %>% select(-Name, -Description) %>% as.matrix
rownames(mat) <- df$Name

#
# run xCell
#
res <- xCellAnalysis(mat, parallel.sz=threads)

#
# extract and save cell fractions
#
cf <- data.frame(Name=rownames(res), Description="xCell", res,
  row.names=NULL, check.names=F)
cf %>% write_gct(outfile)

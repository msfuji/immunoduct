source("scripts/common.R")
library(EPIC)


args <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]

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
# run EPIC
#
res <- EPIC(mat)

#
# extract and save cell fractions
#
cf <- res$cellFractions %>% t
cf <- data.frame(Name=rownames(cf), Description="EPIC", cf,
  row.names=NULL, check.names=F)
cf %>% write_gct(outfile)

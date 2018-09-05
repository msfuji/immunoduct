source("scripts/common.R")
source("scripts/CIBERSORT.R")

args <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]

df <- read_gct(infile)
tmpfile <- tempfile()

#
# save as a matrix
#
df %>% select(-Description) %>% rename(GeneSymbol=Name) %>% write_tsv(tmpfile)

#
# run CIBERSORT
#
res <- CIBERSORT("scripts/LM22.txt", tmpfile, perm=0, QN=F, absolute=F)

#
# extract and save cell fractions
#
cf <- data.frame(Name=rownames(res), Description="xCell", res,
  row.names=NULL, check.names=F)
cf %>% write_gct(outfile)

source("scripts/common.R")
library(MCPcounter)


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
# run MCPcounter
#
genes <- read.table("data/MCPcounter_genes.txt", sep="\t", header=TRUE,
  stringsAsFactors=FALSE, colClasses="character", check.names=FALSE)
res <- MCPcounter.estimate(mat, featuresType="HUGO_symbols", genes=genes)

#
# extract and save cell fractions
#
cf <- data.frame(Name=rownames(res), Description="MCPcounter", res,
  row.names=NULL, check.names=F)
cf %>% write_gct(outfile)

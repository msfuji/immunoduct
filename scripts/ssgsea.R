# Run ssGSEA.
# Masashi Fujita, Aug. 23, 2018

source("scripts/common.R")
library(GSVA)
library(GSEABase)

#
# parse arguments
#
args <- commandArgs(trailingOnly=T)
if(length(args) != 4) {
  stop("[ERROR] invalid number of arguments")
}
threads <- as.integer(args[1])
expr_file <- args[2]
gmt_file <- args[3]
outfile <- args[4]

#
# load and format expression
#
df <- read_gct(expr_file)
uniq_expr <- df %>% group_by(Name) %>% dplyr::slice(1) %>% ungroup
mat <- uniq_expr %>% dplyr::select(-Name, -Description) %>% as.matrix
rownames(mat) <- uniq_expr$Name

#
# load gene sets
#
gsc <- getGmt(gmt_file)

#
# run ssGSEA
#
sig <- gsva(mat, gsc, method="ssgsea", tau=0.75, ssgsea.norm=F,
  parallel.sz=threads)
sig <- data.frame(Name=rownames(sig), sig, check.names=F)

#
# add description of gene sets
#
gs_names <- data.frame(Name=lapply(gsc, setName) %>% unlist,
  Description=lapply(gsc, description) %>% unlist)
df <- gs_names %>% inner_join(sig)

#
# save
#
df %>% write_gct(outfile)

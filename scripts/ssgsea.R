# Run ssGSEA.
# Masashi Fujita, Aug. 23, 2018

source("scripts/common.R")
library(GSVA)
library(GSEABase)

#
# parse arguments
#
if(F){
args <- commandArgs(trailingOnly=T)
if(length(args) != 4) {
  stop("[ERROR] invalid number of arguments")
}
threads <- as.integer(args[1])
expr_file <- args[2]
gmt_file <- args[3]
outfile <- args[4]
}
threads <- 1
expr_file <- "example/LUAD_n10.gct"
gmt_file <- "data/Bindea2013.gmt"
outfile <- "hogehoge.gct"

#
# load and format expression
#
df <- read_gct(expr_file)
uniq_expr <- df %>% as_tibble %>% group_by(Name) %>% slice(1) #%>% ungroup
#mat <- uniq_expr %>% select(-Name, -Description) %>% as.matrix

#
# load gene sets
#
#gene_sets <- getGmt(gmt_file)

#
#
#
#df %>% write_gct(outfile)

# Extract specified rows from GCT file.
# Masashi Fujita, Aug. 22, 2018

source("scripts/common.R")

#
# parse arguments
#
args <- commandArgs(trailingOnly=T)
if(length(args) != 3) {
  stop("[ERROR] invalid number of arguments")
}
expr_file <- args[1]
name_file <- args[2]
outfile <- args[3]

#
# load files
#
expr <- read_gct(expr_file)
name <- fread(name_file, header=F, col.names="Name")

#
# filter
#
df <- expr %>% inner_join(name)

#
# save
#
df %>% write_gct(outfile)

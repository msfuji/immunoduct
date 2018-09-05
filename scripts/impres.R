source("scripts/common.R")


#args <- commandArgs(trailingOnly=T)
infile <- "example/LUAD_n10.gct" #args[1]
outfile <- "hoge.gct" #args[2]

df <- read_gct(infile)

#
# load IMPRES genes
#
impfile <- "data/IMPRES_genes.txt"
imp <- fread(impfile)
imp_genes <- c(imp$Gene1, imp$Gene2) %>% unique

#
# Extract expression of IMPRES genes
#
df <- df %>% filter(Name %in% imp_genes)

#
# remove redundant genes
#
df <- df %>% group_by(Name) %>% slice(1) %>% ungroup

#
# check IMPRES genes
#
missing <- imp_genes[!(imp_genes %in% df$Name)]
if(length(missing) > 0){
  print(missing)
  stop("[ERROR] some IMPRES genes were not found in the input.")
}

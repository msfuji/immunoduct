source("scripts/common.R")
library(tidyr)

args <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]

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
  stop("[ERROR] signature gene was not found in the input.")
}

df <- df %>% select(-Description) %>%
  gather(key=sample, value=expression, -Name)

imp.df <- imp %>%
  inner_join(df, by=c("Gene1"="Name")) %>%
  inner_join(df, by=c("Gene2"="Name", "sample"), suffix=c("1", "2"))

#
# compute IMPRES
#
res <- imp.df %>% group_by(sample) %>%
  summarize(impres=sum(expression1 > expression2))

#
# save
#
v <- res$impres
names(v) <- res$sample
t <- data.frame(Name="IMPRES", Description="Signature", t(v), check.names=F)
t %>% write_gct(outfile)

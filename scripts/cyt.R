source("scripts/common.R")


args <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]

df <- read_gct(infile)

human_genes <- c("GZMA", "PRF1")
mouse_genes <- c("Gzma", "Prf1")

if(all(human_genes %in% df$Name)) {
  df <- df %>% filter(Name %in% human_genes)
} else if(all(mouse_genes %in% df$Name)) {
  df <- df %>% filter(Name %in% mouse_genes)
} else {
  msg <- paste("[ERROR] Either GZMA or PRF1 was not found in", infile)
  stop(msg)
}

#
# remove redundant genes
#
df <- df %>% group_by(Name) %>% slice(1) %>% ungroup

df <- df %>% select(-Name, -Description)

cyt <- sqrt(df[1,] * df[2,])
cyt <- data.frame(
  Name="Cytolytic activity",
  Description="Signature",
  cyt, check.names=F
)

cyt %>% write_gct(outfile)

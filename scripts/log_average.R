source("scripts/common.R")
library(preprocessCore)

# do quantile normalization?
qn <- F

# pseudocount for log transformation
offset <- 0.01

args <- commandArgs(trailingOnly=T)
if(args[1]=="--qn") {
  qn <- T
  args <- args[-1]
}
expr_file <- args[1]
gmt_file <- args[2]
outfile <- args[3]

#
# load files
#
df <- read_gct(expr_file)
gsc <- getGmt(gmt_file)

#
# remove redundant genes
#
df <- df %>% group_by(Name) %>% slice(1) %>% ungroup

#
# make matrix
#
mat <- df %>% select(-Name, -Description) %>% as.matrix
rownames(mat) <- df$Name

#
# quantile normalization
#
if(qn) {
  cols <- colnames(mat)
  rows <- rownames(mat)
  mat <- normalize.quantiles(mat)
  colnames(mat) <- cols
  rownames(mat) <- rows
}

#
# log transformation
#
mat <- log10(mat + offset)

#
# average expression
#
res <- lapply(gsc, function(gs) {
  sig_genes <- geneIds(gs)
  sig_mat <- mat[sig_genes,]
  missing <- sig_genes[!(sig_genes %in% rownames(sig_mat))]
  if(length(missing) > 0){
    print(missing)
    stop("[ERROR] signature gene was not found in the input.")
  }
  mean_vec <- apply(sig_mat, 2, mean)
  row <- data.frame(Name=setName(gs), Description=description(gs),
    t(mean_vec), check.names=F, stringsAsFactors=F)
  return(row)
})
res <- do.call(rbind, res)

#
# save
#
res %>% write_gct(outfile)

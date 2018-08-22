library(dplyr)
library(data.table)
library(readr)

#
# function for reading GCT
#
read_gct <- function(filename) {
  df <- fread(filename, skip=2, sep="\t")

  # check file format
  cols <- colnames(df)
  if(cols[1] != "Name" || cols[2] != "Description") {
    msg <- sprintf("[ERROR] invalid column names: %s", filename)
    stop(msg)
  }
  return(df)
}

#
# function for writing GCT
#
write_gct <- function(df, outfile) {
  "#1.2\n" %>% cat(file=outfile)
  sprintf("%d\t%d\n", nrow(df), ncol(df)-2) %>% cat(file=outfile, append=T)
  df %>% write_tsv(outfile, append=T, col_names=T)
}

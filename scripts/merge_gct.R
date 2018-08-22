# Merge multiple GCT files into one.
# Masashi Fujita, Aug. 22, 2018

library(dplyr)
library(data.table)
library(readr)

#
# parse arguments
#
args <- commandArgs(trailingOnly=T)
if(length(args) < 3) {
  stop("[ERROR] too few input arguments")
}
dir <- args[1]  # direction to merge tables
outfile <- args[2]
infiles <- args[c(-1, -2)]
if(dir != "row" && dir != "col") {
  stop("[ERROR] The 1st argument must be either 'row' or 'col'.")
}

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

#
# To check consistency of Name and Description between two GCTs
#
check_names <- function(df1, df2) {
  same_name <- all(df1$Name == df2$Name)
  same_desc <- all(df1$Description == df2$Description)
  return (same_name && same_desc)
}

#
# combine GCTs
#
if(length(infiles) == 1) {
  file.copy(infiles[1], outfile)
} else {
  # bind columns
  if(dir == "col") {
    dfs <- infiles %>% lapply(read_gct)
    for(i in 2:length(dfs)) {
      if(!check_names(dfs[[1]], dfs[[i]])) {
        stop("[ERROR] Incosistent Name or Description between GCTs")
      }
      dfs[[i]] <- dfs[[i]] %>% select(-Name, -Description)
    }
    df <- bind_cols(dfs)
  }
  # bind rows
  else {
    df <- infiles %>% lapply(read_gct) %>% bind_rows
  }
  df %>% write_gct(outfile)
}

# Merge multiple GCT files into one.
# Masashi Fujita, Aug. 22, 2018

source("scripts/common.R")

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

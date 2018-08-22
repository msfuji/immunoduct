# Combine GCT files and FPKM-UQ from raw read counts.
# Masashi Fujita, 1208/2017

library(dplyr)
library(data.table)
library(readr)

args <- commandArgs(trailingOnly=T)
count_file <- args[1]
name_file <- args[2]
outdir <- args[3]

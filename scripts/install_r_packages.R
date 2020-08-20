install.packages(c("dplyr", "tidyr", "data.table", "readr", "RCurl", "bit64", "e1071", "lazyeval", "prettyunits", "backports", "markdown"))
install.packages("https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz", repos=NULL, type="source")

install.packages("devtools")

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA")
BiocManager::install("preprocessCore")

# ESTIMATE
install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)

# EPIC
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

# MCPcounter
devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")

# xCEll
options(unzip = "internal")
devtools::install_github('dviraran/xCell')

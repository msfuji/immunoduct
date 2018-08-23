# ESTIMATE
install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)

# EPIC
install.packages("devtools")
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

# MCPcounter
install.packages("devtools","curl") ##Installs devtools and the MCPcounter dependancy 'curl'
library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")

# xCEll
devtools::install_github('dviraran/xCell')

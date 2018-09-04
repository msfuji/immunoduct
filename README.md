# immunoduct
Pipeline for tumor immunology. Input should be GCT files made by `riboduct`.

## Installation
### 1. Install conda
`immunoduct` requires `conda` package manager. To install `conda` for Linux,
```
wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
bash Anaconda3-5.2.0-Linux-x86_64.sh
echo ". ${HOME}/anaconda3/etc/profile.d/conda.sh" >> ~/.bashrc
source ~/.bashrc
```

### 2. Install pipeline
Download and install `immunoduct`.
```
git clone https://github.com/msfuji/immunoduct.git
cd immunoduct
conda activate
bash immunoduct.sh install
conda deactivate
```

### 3. Install additional software
Install `clustergrammer` for visualizatoin. Also install `ESTIMATE`, `xCell`,
`EPIC`, and `MCPcounter` for immune profiling.
```
conda activate immunoduct

# install clustergrammer
pip install clustergrammer_widget
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter nbextension enable --py --sys-prefix clustergrammer_widget

# install ESTIMATE, MCPcounter, xCELL, EPIC
Rscript scripts/install_r_packages.R

conda deactivate
```

`CIBERSORT` is not automatically installed because of its license.
While `immunoduct` runs without `CIBERSORT`, you may go to
https://cibersort.stanford.edu/ and request R source code.
Download `CIBERSORT.R`, `LM22.txt`, and place them in `scripts/`.

## Usage
Make a local copy of pipeline for each project.
```
git clone https://github.com/msfuji/immunoduct.git
cd immunoduct
```
Modify `config.yaml`. Start running the pipeline.
```
conda activate immunoduct
bash immunoduct.sh run config.yaml
conda deactivate
```
The output is saved to `output/immunoduct.gct`.

## Visualization
```
conda activate immunoduct
bash immunoduct.sh view
conda deactivate
```

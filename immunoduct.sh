#!/bin/bash

command=$1
config_yaml=$2

#
# install immunoduct
#
if [ $command = "install" ]; then
  echo "Installing immunoduct..."
  conda config --add channels r
  conda config --add channels defaults
  conda config --add channels bioconda
  conda env create --name immunoduct --file environment.yaml
  echo "DONE."
  exit 0
elif [ $command = "view" ]; then
  jupyter notebook cluster/cluster.ipynb
elif [ -z $config_yaml ]; then
  echo "[ERROR] config.yaml not specified."
  exit 1
fi

#
# check config.yaml whether to use SGE
#
sge_jobs=`awk '$1~/^sge_jobs:/ {print $2}' config.yaml`
if [ -z $sge_jobs ]; then
  use_sge=false
  echo "Use SGE cluster... [No]"
else
  use_sge=true
  echo "Use SGE cluster... [Yes]"
  echo "Max SGE jobs... ["$sge_jobs"]"
fi

#
# add SGE options for snakemake
#
snake_command="snakemake --config env_dir=$CONDA_PREFIX --configfile $config_yaml"
if $use_sge; then
  snake_command=`echo $snake_command --cluster \"qsub -pe def_slot {threads} -o {log} -e {log}\" --jobs $sge_jobs`
fi

#
# run snakemake
#
if [ $command = "run" ]; then
  echo "Running immunoduct..."
  bash -c "$snake_command"
  echo "DONE."
else
  echo "[ERROR] Unknown command ", $command
fi

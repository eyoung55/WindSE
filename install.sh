#!/bin/bash 

### Run this to use conda in the script
# eval "$(conda shell.bash hook)"
source $(conda info --base)/etc/profile.d/conda.sh

### Create the Environment
conda create -y --name $1
conda activate $1

### Install conda-forge dependencies 
conda install -y -c conda-forge fenics=2019.1.0 dolfin-adjoint matplotlib scipy mshr hdf5 pyyaml memory_profiler pytest

### Install the new tsfc compiler
pip install git+https://github.com/blechta/tsfc.git@2018.1.0
pip install git+https://github.com/blechta/COFFEE.git@2018.1.0
pip install git+https://github.com/blechta/FInAT.git@2018.1.0
pip install singledispatch networkx pulp

### Install editible version of WindSE
pip install -e .

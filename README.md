# Description

## Installation

Prerequisites to installation :
- miniconda (https://docs.anaconda.com/free/miniconda/miniconda-install/) 
- git (https://git-scm.com/downloads)

First, in a terminal, clone this repository and its unpackaged dependancies in the desired location :
```
git clone https://github.com/GeraultTr/Root_BRIDGES.git
git clone https://github.com/GeraultTr/genericmodel.git
git clone https://github.com/GeraultTr/data_utility.git
git clone https://github.com/GeraultTr/Root_CyNAPS.git
git clone https://github.com/GeraultTr/Rhizodep.git
cd root_bridges
```

Then, install packaged dependancies with : 
```
conda env create -f requirements.yml
conda activate root_bridges
```

Finally, setup the locally unpackaged dependancies : 
```
cd .../root_bridges
python -m setup.py develop
cd .../rhizodep
python -m setup.py develop
cd .../Root_CyNAPS
python -m setup.py develop
cd .../genericmodel
python -m setup.py develop
cd .../data_utility
python -m setup.py develop
```
# Electronic Notebook: Generalized Moment Correction for Long-Ranged Electrostatics

This repository contains detailed documentation for reproducing the data and analysis presented in the research paper

**Generalized Moment Correction for Long-Ranged Electrostatics**  
Björn Stenqvist, Vidar Aspelin, and Mikael Lund  
_Journal of Chemical Theory and Computation_, 2020  
DOI: [10.1021/acs.jctc.9b01003](https://doi.org/10.1021/acs.jctc.9b01003)  
 
## Layout

- `bulk/` Molecular Dynamics setup
- `salts/` Simulations of aqueous solutions of NaCl and NaI, Kirkwood-Buff analysis

## Usage

To open the Notebooks, install python via [Miniconda](https://conda.io/miniconda.html) and
make sure all required packages are loaded by issuing the following terminal commands

``` bash
    conda env create -f environment.yml
    source activate qpotential-*
    jupyter-notebook
```

where `*` can be either `bulk` or `salts` depending on the subdirectory above.

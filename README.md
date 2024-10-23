[![DOI](https://zenodo.org/badge/688868849.svg)](https://zenodo.org/badge/latestdoi/688868849)

This repository contains code that reproduces the main analyses of the following study:

>Dimitrios - Georgios Kontopoulos, Arnaud Sentis, Martin Daufresne, Natalia Glazman, Anthony I. Dell, and Samraat Pawar: <a href='https://doi.org/10.1038/s41467-024-53046-2'>**No universal mathematical model for thermal performance curves across traits and taxonomic groups**</a>. <i>Nature Communications</i> 15, 8855 (2024).

---
 
#### Execution

Each script can be run from the command line or from within R using the 
source command.

The very first script is "TPC_fitting_pipeline.R". It has to be run from the command line as follows:
```
./TPC_fitting_pipeline.R traits # models fitted to trait datasets
```
or
```
./TPC_fitting_pipeline.R enzymes # models fitted to enzyme datasets
```

For this script to run, the [data](https://doi.org/10.6084/m9.figshare.24106161.v2)
need to be placed under a "Data/" directory, outside "Code/". The "Results/model_fits/" 
and "Results/enzyme_model_fits/" directories should also exist.

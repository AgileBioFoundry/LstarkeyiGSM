# LstarkeyiGSM (iLst996)
Lipomyces starkeyi code for developing GSM iLst996 based on strain NRRL Y-11557 genome and NRRL Y-11558 phenotypic data. 

- [Setup](#setup)
- [Summary of files](#summary-of-files)
- [Reference](#reference)
- [License](#license)


## Setup

Model development utilizes the cobra toolbox to perform simulations and has the following dependences:

* [cobra](https://opencobra.github.io/cobrapy/)
* scipy
* pandas
* numpy
* matplotlib
* [scikit-learn](https://scikit-learn.org/stable/index.html)
* seaborn

If the required packages are not installed, you can use the following line from the command prompt to install:
> pip install cobra scipy pandas numpy matplotlib seaborn scikit-learn

## Summary of files

* BiologData_experimental - [Biolog phenotypic plates](https://www.biolog.com/) experimental data and analysis files. 
* Notebooks - notebooks used to build the GSM, with seven seperate steps (labeled step1-step7). Also includes notebooks for simulating GSM fluxes and for essential gene examinations. 
* blastp - results from using ncbi protein blast of the _L. starkeyi_ NRRL Y-11557 and _R. sporidium_ IFO0880_4. 
* data - Omic data, lipidData, biolog summary data, and Lipst1_1genome used for building and refining the GSM.
* finalGSM - iLst996 in various formats.
* memoteReport - memoteReport for iLst996.
* models - GSMs used and generated in this work. 
* orthoMCL - load and read the orthoMCL results performed on _L. starkeyi_ NRRL Y-11557, _R. sporidium_ IFO0880_4, _S. cerevisiae_, and _Y. lipolytica_.
* orthoMCL_multispecies - orthoMCL results and analysis for _Lipomyces_ clade and outlying species (26 species). 

## Reference

This project has been published in **Genome-scale model development and genomic sequencing of the oleaginous clade Lipomyces**. _Frontiers in Bioenergy and Biotechnology_, Volume 12. DOI: 10.3389/fbioe.2024.1356551. Link to [online article](https://www.frontiersin.org/articles/10.3389/fbioe.2024.1356551/full).



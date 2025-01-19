# EBG: Educated Bootstrap Guesser
[![image](https://img.shields.io/pypi/v/ebg.svg)](https://pypi.python.org/pypi/ebg)
[![image](https://img.shields.io/conda/vn/conda-forge/ebg.svg)](https://anaconda.org/conda-forge/ebg)
[![image](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/GPL-3-0)

Documentation: https://github.com/wiegertj/EBG/wiki
## Description

EBG is a Python tool for predicting the Felsenstein Bootstrap Support of phylogenies inferred by RAxML-NG.
It was trained on empirical datasets from TreeBASE and can use both AA and DNA data.


## Installation
### Using conda
The latest version of EBG can easily be installed via conda:
```
conda install ebg -c conda-forge
```
### Using pip
```
pip install ebg
```
## Usage Example - Light Mode
EBG now offers a faster, lightweight alternative to the standard EBG mode by omitting the most computationally expensive feature calculations while maintaining performance comparable to the standard mode. This mode is particularly useful for MSAs with a large number of sites.

In our experiments, performance decreases by 6% in light mode.

You can enable light mode in EBG using the ```-light``` flag.

## Usage Example - Standard Mode
A simple command line call of EBG looks like this:
```
ebg -msa /test/example.fasta -tree /test/example.bestTree -model /test/example.bestModel -t b -o test 
```
This command will use the MSA in fasta format, and the best tree inferred with RAxML-NG and the model.
By selecting ```-t b```(oth) EBG will output the bootstrap predictions as well as the probabilities for exceeding different bootstrap thresholds (70, 75, 80, 85). 
The results will be stored in a folder called test.

Please keep in mind that EBG requires an installation of RAxML-NG. By default, it uses the command ```raxml-ng```. 
If your RAxML-NG installation is not part of the PATH variable, you can specify the path to the RAxML-NG binary file with the parameter ```-raxmlng PATH_TO_RAXMLNG```.

### Citation
If you are using EBG for your publication, please cite our published paper in Molecular Biology and Evolution: 
[EBG Paper](https://academic.oup.com/mbe/article/41/10/msae215/7825466)
### References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
[https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

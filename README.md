# EBG: Educated Bootstrap Guesser

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
## Usage Example
A simple command line call of EBG looks like this:
```
ebg -msa /test/example.fasta -tree /test/example.bestTree -model /test/example.bestModel -t b -o test 
```
This command will use the MSA in fasta format, and the best tree inferred with RAxML-NG and the model.
By selecting ```-t b```(oth) EBG will output the bootstrap predictions as well as the probabilities for exceeding different bootstrap thresholds (70, 75, 80, 85). 
The results will be stored in a folder called test.

Please keep in mind that EBG requires an installation of RAxML-NG. By default, it uses the command ```raxml-ng```. 
If your RAxML-NG installation is not part of the PATH variable, you can specify the path to the RAxML-NG binary file with the parameter ```-raxmlng PATH_TO_RAXMLNG```.
### References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
[https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

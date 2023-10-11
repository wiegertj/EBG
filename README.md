# EBG: Educated Bootstrap Guesser

## Description

EBG is a python tool for predicting the Felsenstein Bootstrap Support of phylogenies inferred by RAxML-NG.
It was trained on empirical datasets from the TreeBASE and is able to use both AA and DNA data.


## Installation
### Using conda
The latest version of EBG can easily be installed via conda:
```
conda install ebg -c bioconda
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
By selecting ```-t b```(oth) EBG will ouput the bootstrap predictions as well as the probabilities for exceeding different bootstrap thresholds (70, 75, 80, 85). 
The results will be stored in a folder called test.
### References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
[https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

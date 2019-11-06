# Sphetcher

Sphetcher is a C++ library to downsample singe-cell RNA secquencing datasets. 


# Installation
All dependencies are bundled with Sphetcher. To build Sphetcher, go to directory src/, and simply run

```
make
```

## Running Sphetcher ##

### To begin ###

First you will need to provide the gene expression matrix in separate comma format (.csv) where each rows are samples and columns are features (genes, transcripts).

additionally you can provide the prior information (e.g, labels, time of the experiments) in case of fair sampling:

An example of inputs is provided in the directory ```/data```. 

### Usage ###

Once you have compiled Sphetcher it can be run easily with the following command:

```
sphetcher expression_matrix.csv
```

#### Input/Output formats

Input: 

`expression_matrix.csv`
  : expression matrix: rows are samples, columns are genes.

Output:

Indices of the sketch: __INDEX START FROM ZERO (0) __



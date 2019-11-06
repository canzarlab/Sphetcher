# Sphetcher

Sphetcher is a C++ library to sample big singe-cell RNA secquencing datasets based on the spherical thresholding algorithm. It selects a small subset of cells referred to as ```sketch``` that evenly cover the transcriptomic space occupied by the original dataset. Such a sketch can accelerate downstream analyses and highlight rare cell types.


# Installation
All dependencies are bundled with Sphetcher. To build Sphetcher, go to directory src/, and run

```
make
```

## Running Sphetcher ##

### To begin ###

First you will need to provide the gene expression matrix in separate comma format (.csv), where rows are samples and columns are features (genes, transcripts).

Additionally you can provide the prior information (e.g, cell label, collection time point) in case you want to preserve certain number of samples from each category. 

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

Indices of the sketch: __INDICES START FROM ZERO (0)__



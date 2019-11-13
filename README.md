# Sphetcher

A software package for sampling massive singe-cell RNA secquencing datasets based on the spherical thresholding algorithm. It selects a small subset of cells referred to as ```sketch``` that evenly cover the transcriptomic space occupied by the original dataset. Such a sketch can accelerate downstream analyses and highlight rare cell types.


# Installation
A compiler that supports C++11 is needed to build sphetcher. You can download and compile the latest code from github as follows:

```
git clone  https://github.com/canzarlab/Sphetcher
cd src
make
```

## Running Sphetcher ##

### To begin ###

First you will need to provide the gene expression matrix in separate comma format (.csv), where rows are samples (cells) and columns are features (genes, transcripts, principal components PCs).

Additionally you can provide the prior information (e.g, cell label, collection time point) in case you want to preserve certain number of samples from each category (`l_min`). 

An example of inputs is provided in the directory ```/data```. 

### Usage ###

Once you have compiled Sphetcher it can be run easily with one of the following two options:

```
sphetcher expression_matrix.csv sketch_size
```
or 
```
sphetcher expression_matrix.csv sketch_size class_labels.csv l_min
```
For an example provided in ```/data```
```
sphetcher zeisel_pca.csv 1000 
or 
sphetcher zeisel_pca.csv 1000 zeisel_pca_labels.csv 3 
```

#### Input/Output formats

Input: 

`expression_matrix.csv`
  : expression matrix: rows are samples, columns are genes. <br/>
 `sketch_size` 
  : number of samples to obtain from the data set. <br/>
`class_labels.csv`
  : prior information, each class is presented by an integer between 1 and K, where K is # of the classes. <br/>
`l_min`
  : minimum number of samples we want to sample from each class <br/>

Output:

An indicator vector (`indicator_solutions.csv`) of size ```n```(total # of samples) where 1 indicates the sample is in the sketch, 0 otherwise. 



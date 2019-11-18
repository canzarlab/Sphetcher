# Sphetcher
A software package for sampling massive singe-cell RNA secquencing datasets based on the spherical thresholding algorithm. It selects a small subset of cells referred to as ```sketch``` that evenly cover the transcriptomic space occupied by the original dataset. Such a sketch can accelerate downstream analyses and highlight rare cell types.

<img src=img/overviewv2.png  width="100%" height = "550">

# Installation
A compiler that supports C++11 is needed to build sphetcher. You can download and compile the latest code from github as follows:

```
git clone  https://github.com/canzarlab/Sphetcher
cd src
make
```

## Running Sphetcher ##

### To begin ###

First you need to provide a matrix where rows are samples (cells) and columns are features (genes, transcripts, principal components PCs).

Additionally you can provide the prior information (e.g, cell label, collection time point) in case you want to preserve certain number of samples from each category. 

An example of inputs is provided in the directory ```/data```. 

### Usage ###

Once you have compiled Sphetcher it can be run easily with one of the following two options:

```
sphetcher expression_matrix.csv sketch_size sketch_indicator_output.csv
```
or 
```
sphetcher expression_matrix.csv sketch_size class_labels.csv l_min sketch_indicator_output.csv
```
For an example provided in ```/data```
```
sphetcher zeisel_pca.csv 1000 sketch_indicator_output.csv
or 
sphetcher zeisel_pca.csv 1000 zeisel_pca_labels.csv 3 sketch_indicator_output.csv
```

#### Input/Output formats

Input: 

`expression_matrix.csv`
  : expression matrix in comma-separated values (CSV) format: rows are cells, columns are features. <br/>
 `sketch_size` 
  : number of samples to obtain from the data set. <br/>
`class_labels.csv`
  : prior information stored in a column vector, each class is presented by an integer between 1 and K, where K is the number of classes. <br/>
`l_min`
  : minimum number of representatives we want to sample from each class <br/>

Output:

`sketch_indicator_output.csv` : an indicator vector of `n` samples where 1 indicates the sample is in the sketch, 0 otherwise (there are `sketch_size` 1's in the vector).


